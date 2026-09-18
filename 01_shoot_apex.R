# Submission edition: run from this directory or set POPLAR_CODE_ROOT.
.code_root <- Sys.getenv("POPLAR_CODE_ROOT", unset = "")
if (!nzchar(.code_root)) {
  .script <- grep("^--file=", commandArgs(), value = TRUE)
  .code_root <- if (length(.script)) dirname(normalizePath(sub("^--file=", "", .script[1]))) else getwd()
}
source(file.path(.code_root, "R", "submission_setup.R"))

stages <- list(
  preprocess = function() {

source(file.path(ATLAS_CODE_ROOT, "R", "common_helpers.R"), local = TRUE)
init_project_dirs()

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(clusterProfiler)
  library(RcppML)
  library(qs)
  library(scCustomize)
  library(harmony)
})

make_output_dirs()
plot_cfg <- list(resolution_width = 14, resolution_height = 12, final_width = 16, final_height = 10, marker_width = 16, marker_height = 12, ncol = 4, pt.size.factor = 1.8, crop = FALSE)

# --------------------------------------- #
# SAM pipeline (section-splitting version)
# --------------------------------------- #

# --------------------------------------- #
# 0. Manually remove damaged SAM regions ##
# --------------------------------------- #

sam_cleanup_images <- c("sam_B_s1", "sam_B_s2", "sam_C_s3")
sam_cleanup_polygon_file <- "tables/SAM_detached_leaf_exclusion_polygons.csv"
sam_cleanup_removed_file <- "tables/SAM_manual_cleanup_removed_spots.csv"
redraw_sam_cleanup_polygons <- FALSE

get_sam_section_coordinates <- function(obj, image) {
  if (!image %in% names(obj@images)) stop("Image not found: ", image)
  coords <- as.data.frame(GetTissueCoordinates(obj, image = image, scale = NULL))
  x_col <- intersect(c("x", "imagecol", "col"), names(coords))[1]
  y_col <- intersect(c("y", "imagerow", "row"), names(coords))[1]
  if (is.na(x_col) || is.na(y_col)) stop("Could not identify x/y columns for ", image)

  tibble(
    spot = rownames(coords),
    image = image,
    x = coords[[x_col]],
    y = coords[[y_col]]
  ) %>%
    filter(spot %in% colnames(obj), is.finite(x), is.finite(y))
}

point_in_polygon_base <- function(px, py, poly_x, poly_y) {
  inside <- rep(FALSE, length(px))
  j <- length(poly_x)
  for (i in seq_along(poly_x)) {
    crosses <- ((poly_y[i] > py) != (poly_y[j] > py)) &
      (px < (poly_x[j] - poly_x[i]) * (py - poly_y[i]) /
         (poly_y[j] - poly_y[i] + .Machine$double.eps) + poly_x[i])
    inside <- xor(inside, crosses)
    j <- i
  }
  inside
}

select_sam_exclusion_polygon_interactive <- function(
    obj,
    image,
    polygon_id = paste0("cleanup_", image),
    polygon_file = "tables/SAM_detached_leaf_exclusion_polygons.csv",
    append = FALSE
) {
  if (!interactive()) stop("Run this function in an interactive RStudio session.")
  coords <- get_sam_section_coordinates(obj, image)

  plot(
    coords$x, coords$y,
    ylim = rev(range(coords$y)),
    asp = 1,
    pch = 1,
    cex = 0.45,
    col = "grey40",
    xlab = "Spatial x",
    ylab = "Spatial y",
    main = paste0(
      image, " (", nrow(coords), " spots)\n",
      "Click around region to remove; press Esc/Finish when done"
    )
  )

  selected <- locator(type = "l", col = "#D73027", lwd = 2)
  if (is.null(selected) || length(selected$x) < 3) {
    stop("At least three polygon vertices are required.")
  }

  polygon <- tibble(
    image = image,
    polygon_id = polygon_id,
    vertex_order = seq_along(selected$x),
    x = selected$x,
    y = selected$y
  )
  inside <- point_in_polygon_base(coords$x, coords$y, polygon$x, polygon$y)
  points(coords$x[inside], coords$y[inside], pch = 16, cex = 0.7, col = "#D73027")
  polygon_closed <- bind_rows(polygon, polygon[1, ])
  lines(polygon_closed$x, polygon_closed$y, col = "#D73027", lwd = 2)
  title(sub = paste(sum(inside), "spots selected"))

  if (append && file.exists(polygon_file)) {
    existing <- read.csv(polygon_file, check.names = FALSE)
    polygon <- bind_rows(existing, polygon)
  }
  write.csv(polygon, polygon_file, row.names = FALSE)
  invisible(polygon)
}

select_sam_cleanup_polygons_interactive <- function(
    obj,
    images = c("sam_B_s1", "sam_B_s2", "sam_C_s3"),
    polygon_file = "tables/SAM_detached_leaf_exclusion_polygons.csv"
) {
  missing_images <- setdiff(images, names(obj@images))
  if (length(missing_images) > 0) {
    stop("Images not found: ", paste(missing_images, collapse = ", "))
  }

  for (i in seq_along(images)) {
    message("Selecting cleanup region ", i, "/", length(images), ": ", images[[i]])
    select_sam_exclusion_polygon_interactive(
      obj = obj,
      image = images[[i]],
      polygon_id = paste0("cleanup_", images[[i]]),
      polygon_file = polygon_file,
      append = i > 1
    )
  }
  invisible(read.csv(polygon_file, check.names = FALSE))
}

remove_sam_manual_cleanup_spots <- function(
    obj,
    polygon_file = "tables/SAM_detached_leaf_exclusion_polygons.csv",
    removed_file = "tables/SAM_manual_cleanup_removed_spots.csv",
    max_removed_fraction = 0.5
) {
  if (!file.exists(polygon_file)) {
    message("No manual cleanup polygon file found; returning the unfiltered object.")
    return(obj)
  }

  polygons <- read.csv(polygon_file, check.names = FALSE) %>%
    arrange(image, polygon_id, vertex_order)
  required <- c("image", "polygon_id", "vertex_order", "x", "y")
  if (!all(required %in% names(polygons))) {
    stop("Polygon file must contain: ", paste(required, collapse = ", "))
  }

  keys <- polygons %>% distinct(image, polygon_id)
  removed <- purrr::map2_dfr(keys$image, keys$polygon_id, function(img, poly_id) {
    polygon <- polygons %>%
      filter(.data$image == .env$img, .data$polygon_id == .env$poly_id)
    coords <- get_sam_section_coordinates(obj, img)
    coords %>%
      mutate(
        polygon_id = poly_id,
        remove = point_in_polygon_base(x, y, polygon$x, polygon$y)
      ) %>%
      filter(remove)
  }) %>%
    distinct(spot, .keep_all = TRUE)

  write.csv(removed, removed_file, row.names = FALSE)
  removal_check <- removed %>%
    dplyr::count(image, name = "n_remove") %>%
    mutate(
      n_section = vapply(
        image,
        function(img) nrow(get_sam_section_coordinates(obj, img)),
        integer(1)
      ),
      fraction_removed = n_remove / n_section
    )
  print(removal_check)
  if (any(removal_check$fraction_removed > max_removed_fraction)) {
    stop("Safety stop: more than 50% of a section was selected. Redraw its polygon.")
  }

  qc_coords <- purrr::map_dfr(unique(polygons$image), ~ get_sam_section_coordinates(obj, .x)) %>%
    mutate(removed = spot %in% removed$spot)
  qc_plot <- ggplot(qc_coords, aes(x, y, color = removed)) +
    geom_point(size = 0.55) +
    scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D73027")) +
    scale_y_reverse() +
    coord_fixed() +
    facet_wrap(~image) +
    theme_classic() +
    labs(color = "Removed")
  ggsave("QC/SAM_manual_cleanup_QC.pdf", qc_plot, width = 10, height = 5)

  message("Removing ", nrow(removed), " manually selected SAM spots.")
  subset(obj, cells = setdiff(colnames(obj), removed$spot))
}

# --------------------------------------- #
# 1. Load cleaned raw merged SAM object ##
# --------------------------------------- #
sam_filtered_raw_file <- "saved_obj/sobj_sam_split_raw_manual_cleanup_filtered.qs"
if (!file.exists(sam_filtered_raw_file)) stop("Filtered raw SAM object not found: ", sam_filtered_raw_file)
sobj_sam <- qread(sam_filtered_raw_file)
if (!"Spatial" %in% Assays(sobj_sam)) stop("The filtered SAM object has no Spatial assay.")
count_layers <- grep("^counts", Layers(sobj_sam[["Spatial"]]), value = TRUE)
if (length(count_layers) == 0) stop("The Spatial assay has no raw count layers.")
if (length(count_layers) > 1) sobj_sam <- JoinLayers(sobj_sam, assay = "Spatial")
if (!"counts" %in% Layers(sobj_sam[["Spatial"]])) stop("Failed to join Spatial count layers.")

# Remove stale downstream results if the object was previously processed.
if ("SCT" %in% Assays(sobj_sam)) sobj_sam[["SCT"]] <- NULL
sobj_sam@reductions <- list()
sobj_sam@graphs <- list()
sobj_sam@neighbors <- list()
DefaultAssay(sobj_sam) <- "Spatial"
message("Loaded ", ncol(sobj_sam), " cleaned SAM spots from ", sam_filtered_raw_file)

# --------------------------------------- #
# 2. Normalize and integrate #############
# --------------------------------------- #
options(future.globals.maxSize = 2 * 1024^3)
future::plan("sequential")
sobj_sam <- process_spatial_integration(sobj_sam, resolution = 0.6, dims = 1:30)
# Optional intermediate checkpoint; comment out to reduce disk use.
# qsave(sobj_sam, "saved_obj/sobj_sam_split_cleaned_harmony_res0.6.qs_2026.7.8")
# sobj_sam <- qread("saved_obj/sobj_sam_split_cleaned_harmony_res0.6.qs_2026.7.8")

# --------------------------------------- #
# 4. QC + resolution sweep ##############
# --------------------------------------- #
plot1 <- VlnPlot_scCustom(sobj_sam, features = c("nCount_Spatial", "nFeature_Spatial"), group.by = "orig.ident", plot_median = TRUE) + NoLegend()
plot2 <- SpatialFeaturePlot(sobj_sam, features = c("nCount_Spatial", "nFeature_Spatial"), crop = FALSE, ncol = 4)
pdf("QC/QC_sam_split_violin.pdf", width = 11, height = 7); print(plot1); dev.off()
pdf("QC/QC_sam_split_spatial.pdf", width = 14, height = 10); print(plot2); dev.off()
write.csv(describeBy(sobj_sam@meta.data, group = sobj_sam@meta.data$orig.ident, mat = TRUE), "QC/QC_stats_sam_split.csv", row.names = FALSE)

resolutions <- c(0.4, 0.5, 0.6, 0.8)
pdf("clustering/sam_resolution_optimization_s4.pdf", width = 14, height = 12)
for (res in resolutions) {
  sobj_tmp <- FindClusters(sobj_sam, resolution = res)
  p <- SpatialDimPlot(sobj_tmp, cells.highlight = CellsByIdentities(sobj_tmp), facet.highlight = TRUE,
                      images = names(sobj_tmp@images)[4], ncol = 6, pt.size.factor = 5) +
    plot_annotation(title = paste0("SAM resolution: ", res))
  print(p)
}
dev.off()

## after visual review, set final resolution
### final resolution: 0.4#################
sobj_sam <- FindClusters(sobj_sam, resolution = 0.4)

## (optimoal) adjust subcluster based on the marker and celltype disctribution. 
sobj_sam_sub8 <- FindSubCluster(sobj_sam, cluster = "8", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "sub8")
Idents(sobj_sam_sub8) <- sobj_sam_sub8$sub8
SpatialDimPlot(sobj_sam_sub8, images = names(sobj_tmp@images)[1], pt.size.factor = 1, group.by = 'sub8')
DimPlot(sobj_sam_sub8,group.by = 'sub8')
LinkedDimPlot(sobj_sam_sub8, group.by = 'sub8', image = names(sobj_tmp@images)[3])
## split cluster 5 (mixed epi and cortex)
sobj_sam_sub5 <- FindSubCluster(sobj_sam, cluster = "5", graph.name = "SCT_snn", resolution = 0.3, subcluster.name = "sub5")
Idents(sobj_sam_sub5) <- sobj_sam_sub5$sub5
SpatialDimPlot(sobj_sam_sub5, images = names(sobj_sam_sub5@images)[3], pt.size.factor = 5, group.by = 'sub5')
DimPlot(sobj_sam_sub5,group.by = 'sub5')
LinkedDimPlot(sobj_sam_sub5, group.by = 'sub5', image = names(sobj_sam_sub5@images)[7])

## put the subcluster back to the main cluster annotation, and rename the subclusters as needed.
sobj_sam <- sobj_sam_sub5
sobj_sam$seurat_clusters <- Idents(sobj_sam)

## after finalized the cluster, plot split cluster again
img_names <- names(sobj_sam@images)[c(1,3,6,8)]  # only plot the 3 split sections, not the shared one
cluster_cells <- CellsByIdentities(sobj_sam)

pdf("clustering/SAM_split_clusters_final_res0.4_2026.7.8.pdf", width = 14, height = 13)
for (img in img_names) {
  p <- SpatialDimPlot(
    sobj_sam,
    cells.highlight = cluster_cells,images = img,
    facet.highlight = TRUE,alpha = 0.8,ncol = 4,
    pt.size.factor = 1.5,crop = FALSE
  ) +
    plot_annotation(title = paste("Final split cluster distribution:", img))
  print(p)
}
dev.off()

# --------------------------------------- #
# 5. Marker-based annotation #############
# --------------------------------------- #
source(file.path(ATLAS_CODE_ROOT, "R", "common_helpers.R"), local = TRUE)
mkr_list <- read.csv("updated_markerlist_poplar_2.8.26.csv")
mkr_list_use <- mkr_list %>% filter(GeneID.v5 %in% rownames(sobj_sam))
make_marker_heatmap(sobj_sam, mkr_list_use, "marker/marker_heatmap_sam_split.pdf",
                    out_table = "marker/marker_heatmap_sam_split_plotted_markers.csv")

## predict the cell type by marker enrichment 
ann_res <- score_and_assign_celltypes(sobj_sam, mkr_list_use, assay = "SCT")
sobj_sam <- ann_res$object
write.csv(ann_res$cluster_scores, "tables/sam_cluster_marker_module_scores.csv", row.names = FALSE)
write.csv(ann_res$top_assign, "tables/sam_cluster_top_predicted_celltypes.csv", row.names = FALSE)
write.csv(ann_res$score_long, "tables/sam_cluster_marker_module_scores_long.csv", row.names = FALSE)

sam_predicted_images <- names(sobj_sam@images)
p_sam_section_umap <- DimPlot(sobj_sam, reduction = "umap", group.by = "orig.ident", alpha = 0.7) + ggtitle("By section")
p_sam_cluster_umap <- DimPlot(sobj_sam, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("By cluster")
p_sam_cluster_spatial <- SpatialDimPlot(sobj_sam, group.by = "seurat_clusters", images = sam_predicted_images, crop = FALSE, ncol = 4, pt.size.factor = 1.6)
p_sam_predicted_spatial <- SpatialDimPlot(sobj_sam, group.by = "predicted_celltype", images = sam_predicted_images, crop = FALSE, ncol = 4, pt.size.factor = 1.6)
pdf("clustering/sam_predicted_annotations_2026.7.8.pdf", width = 16, height = 10)
print(p_sam_section_umap + p_sam_cluster_umap)
print(p_sam_cluster_spatial)
print(p_sam_predicted_spatial)
dev.off()

manual_map_template <- ann_res$top_assign %>%
  transmute(cluster_label, celltype_predicted = predicted_celltype, predicted_score, celltype_final = predicted_celltype) %>%
  arrange(as.numeric(cluster_label))
write_csv(manual_map_template, "tables/manual_cluster_annotation_template_sam.csv")

# --- edit and rerun if needed ---
###  rename manually ######
Idents(sobj_sam) <- sobj_sam$seurat_clusters
sobj_sam_ann <- RenameIdents(sobj_sam, '0' = 'Pith', '1'='Cortex', '2'='Epidermis',
                              '3'='Vasculature','4'='Vasculature', '5_0'='Epidermis','5_1'='Cortex','6'='Leaf primordium', '7'='Leaf primordium','8_0'='Leaf primordium', 
                              '8_1'='Leaf primordium','8_2'='Meristem','9'='Leaf primordium', '10'='Leaf primordium')
sobj_sam_ann$celltypes <- Idents(sobj_sam_ann)
sobj_sam_ann$celltypes_v1 <- as.character(sobj_sam_ann$celltypes)

sam_curated_images <- names(sobj_sam_ann@images)
p_sam_curated_umap_sample <- DimPlot(sobj_sam_ann, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE) + ggtitle("Curated SAM annotation")
p_sam_curated_umap <- DimPlot(sobj_sam_ann, reduction = "umap", group.by = "celltypes", label = TRUE, repel = TRUE) + ggtitle("Curated SAM annotation")
p_sam_curated_spatial <- SpatialDimPlot(sobj_sam_ann, group.by = "celltypes", images = sam_curated_images, crop = FALSE, ncol = 4, pt.size.factor = 1.6)
pdf("clustering/sam_res0.4_annotated_curated_UMAP_2026.7.8.pdf", width = 14, height = 6)
print(p_sam_curated_umap_sample| p_sam_curated_umap)
dev.off()
pdf("clustering/sam_res0.4_annotated_curated_spatial_2026.7.8.pdf", width = 16, height = 9)
print(p_sam_curated_spatial)
dev.off()

### Save annotated final object #######
qsave(sobj_sam_ann, "saved_obj/sobj_sam_split_cleaned_res0.4_annotated_v1_2026.7.17.qs")
# sobj_sam_ann <- qread('saved_obj/sobj_sam_split_cleaned_res0.4_annotated_v1_2026.7.8.qs')

# --------------------------------------- #

sobj_sam_ann <- atlas_annotations(sobj_sam_ann, "sam")
qs::qsave(sobj_sam_ann, "saved_obj/sobj_sam_res0.4_final_v1_2026.7.24.qs")

  },
  markers = function() {
sobj_sam_ann <- qs::qread(atlas_input("saved_obj/sobj_sam_res0.4_final_v1_2026.7.24.qs"))
Seurat::Idents(sobj_sam_ann) <- sobj_sam_ann$celltypes
                                     logfc.threshold = 1, only.pos = TRUE, min.diff.pct = 0.2, min.pct = 0.1)
write.csv(all_de_markers_sam, "marker/sp_all_denovo_markers_sam_split.csv", row.names = FALSE)

top_markers_sam <- Extract_Top_Markers(all_de_markers_sam, num_genes = 5, rank_by = "avg_log2FC",
                                       named_vector = FALSE, make_unique = TRUE, data_frame = TRUE)
pdf("marker/dotplot_sp_top_denovo_markers_sam_split_2026.7.8.pdf", width = 8, height = 10)
print(DotPlot_scCustom(sobj_sam_ann, features = top_markers_sam$gene, flip_axes = TRUE, scale.by = "size",
                       dot.min = 0, dot.scale = 6, facet_label_rotate = TRUE, x_lab_rotate = TRUE))
dev.off()

  },
  rctd = function() {
sobj_sam_ann <- qs::qread(atlas_input("saved_obj/sobj_sam_res0.4_final_v1_2026.7.24.qs"))
Seurat::Idents(sobj_sam_ann) <- sobj_sam_ann$celltypes
# --------------------------------------- #

library(future)
library(qs)
library(dplyr)
library(tibble)
library(tidyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)

plan(sequential)
options(future.globals.maxSize = 30 * 1024^3)

dir.create("saved_obj/scRNA", showWarnings = FALSE, recursive = TRUE)

## Keep Seurat anchor-transfer and RCTD mixture-deconvolution outputs separate.
sam_deconv_dir <- "deconvolution/seurat"
sam_deconv_table_dir <- file.path(sam_deconv_dir, "tables")
sam_deconv_plot_dir <- file.path(sam_deconv_dir, "plots")
sam_rctd_dir <- "deconvolution/RCTD"
sam_rctd_table_dir <- file.path(sam_rctd_dir, "tables")
sam_rctd_plot_dir <- file.path(sam_rctd_dir, "plots")
for (d in c(sam_deconv_dir, sam_deconv_table_dir, sam_deconv_plot_dir, sam_rctd_dir, sam_rctd_table_dir, sam_rctd_plot_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(sam_deconv_dir), dir.exists(sam_rctd_dir))


## 8.1 Load SAM scRNA reference ------------
sobj_scrna_sam <- qread("saved_obj/scRNA/sam_cleaned.qs")

## The SAM scRNA annotation is stored in active identities; copy it to metadata for transfer.
sobj_scrna_sam$scrna_celltype <- as.character(Idents(sobj_scrna_sam))
scrna_celltype_col <- "scrna_celltype"
message("Using SAM scRNA active identities as transfer labels: ", scrna_celltype_col)

table(sobj_scrna_sam[[scrna_celltype_col]][, 1])

scrna_counts_before <- as.data.frame(table(sobj_scrna_sam[[scrna_celltype_col]][, 1]))
colnames(scrna_counts_before) <- c("celltype", "n_cells")
write.csv(scrna_counts_before, file.path(sam_deconv_table_dir, "SAM_scRNA_celltype_counts_before_downsample.csv"), row.names = FALSE)


## 8.2 Plot original SAM scRNA UMAP ------------
sam_scrna_assays <- names(sobj_scrna_sam@assays)
DefaultAssay(sobj_scrna_sam) <- if ("SCT" %in% sam_scrna_assays) "SCT" else if ("RNA" %in% sam_scrna_assays) "RNA" else DefaultAssay(sobj_scrna_sam)

if (!"umap" %in% Reductions(sobj_scrna_sam)) {
  if (!"pca" %in% Reductions(sobj_scrna_sam)) {
    if (!"SCT" %in% names(sobj_scrna_sam@assays)) {
      sobj_scrna_sam <- SCTransform(sobj_scrna_sam, assay = DefaultAssay(sobj_scrna_sam), new.assay.name = "SCT", variable.features.n = 3000, conserve.memory = TRUE, verbose = FALSE)
      DefaultAssay(sobj_scrna_sam) <- "SCT"
    }
    sobj_scrna_sam <- RunPCA(sobj_scrna_sam, npcs = 50, verbose = FALSE)
  }
  sobj_scrna_sam <- RunUMAP(sobj_scrna_sam, reduction = "pca", dims = 1:30, verbose = FALSE)
}

p_sam_scrna_umap <- DimPlot(
  sobj_scrna_sam,
  reduction = "umap",
  group.by = scrna_celltype_col,
  label = TRUE,
  repel = TRUE
) +
  theme_classic() +
  ggtitle("Original SAM scRNA-seq reference")

ggsave(file.path(sam_deconv_plot_dir, "SAM_original_scRNA_UMAP_celltypes.pdf"), p_sam_scrna_umap, width = 9, height = 7)


## 8.3 Proportional downsampling of SAM scRNA reference ------------
set.seed(123)

target_total_cells <- 12000
min_cells_per_celltype <- 50
max_cells_per_celltype <- 1200

scrna_md <- sobj_scrna_sam@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(celltype_use = .data[[scrna_celltype_col]]) %>%
  filter(!is.na(celltype_use))

ct_counts <- scrna_md %>%
  dplyr::count(celltype_use, name = "n_total") %>%
  filter(n_total >= min_cells_per_celltype) %>%
  mutate(
    prop_n = round(n_total / sum(n_total) * target_total_cells),
    n_keep = pmin(max_cells_per_celltype, pmax(min_cells_per_celltype, prop_n), n_total)
  )

print(ct_counts)
write.csv(ct_counts, file.path(sam_deconv_table_dir, "SAM_scRNA_downsample_plan_12k.csv"), row.names = FALSE)

scrna_keep_cells <- scrna_md %>%
  inner_join(ct_counts %>% dplyr::select(celltype_use, n_keep), by = "celltype_use") %>%
  group_by(celltype_use) %>%
  group_modify(~.x[sample(seq_len(nrow(.x)), unique(.x$n_keep)), , drop = FALSE]) %>%
  ungroup() %>%
  pull(cell)

length(scrna_keep_cells)

sobj_scrna_sam_ds <- subset(sobj_scrna_sam, cells = scrna_keep_cells)

scrna_counts_after <- as.data.frame(table(sobj_scrna_sam_ds[[scrna_celltype_col]][, 1]))
colnames(scrna_counts_after) <- c("celltype", "n_cells")
write.csv(scrna_counts_after, file.path(sam_deconv_table_dir, "SAM_scRNA_celltype_counts_after_downsample_12k.csv"), row.names = FALSE)

p_sam_scrna_ds_umap <- DimPlot(
  sobj_scrna_sam_ds,
  reduction = "umap",
  group.by = scrna_celltype_col,
  label = TRUE,
  repel = TRUE
) +
  theme_classic() +
  ggtitle("Downsampled SAM scRNA-seq reference")

ggsave(file.path(sam_deconv_plot_dir, "SAM_downsampled_scRNA_UMAP_celltypes_12k.pdf"), p_sam_scrna_ds_umap, width = 9, height = 7)

# Optional intermediate scRNA reference checkpoint.
# qsave(sobj_scrna_sam_ds, "saved_obj/scRNA/sam_cleaned_downsampled_12k_for_spatial_transfer.qs")


## 8.4 SCT-normalize downsampled SAM scRNA reference ------------
sam_scrna_ds_assays <- names(sobj_scrna_sam_ds@assays)
DefaultAssay(sobj_scrna_sam_ds) <- if ("RNA" %in% sam_scrna_ds_assays) "RNA" else DefaultAssay(sobj_scrna_sam_ds)

if ("SCT" %in% names(sobj_scrna_sam_ds@assays)) sobj_scrna_sam_ds[["SCT"]] <- NULL
gc()

sobj_scrna_sam_ds <- SCTransform(
  sobj_scrna_sam_ds,
  assay = DefaultAssay(sobj_scrna_sam_ds),
  new.assay.name = "SCT",
  variable.features.n = 3000,
  return.only.var.genes = TRUE,
  conserve.memory = TRUE,
  verbose = FALSE
)

DefaultAssay(sobj_scrna_sam_ds) <- "SCT"
sobj_scrna_sam_ds <- RunPCA(sobj_scrna_sam_ds, assay = "SCT", npcs = 50, verbose = FALSE)
sobj_scrna_sam_ds <- RunUMAP(sobj_scrna_sam_ds, reduction = "pca", dims = 1:30, return.model = TRUE, verbose = FALSE)

p_sam_scrna_ds_sct_umap <- DimPlot(
  sobj_scrna_sam_ds,
  reduction = "umap",
  group.by = scrna_celltype_col,
  label = TRUE,
  repel = TRUE
) +
  theme_classic() +
  ggtitle("SCT-normalized downsampled SAM scRNA-seq reference")

ggsave(file.path(sam_deconv_plot_dir, "SAM_downsampled_scRNA_SCT_UMAP_celltypes_12k.pdf"), p_sam_scrna_ds_sct_umap, width = 9, height = 7)

## 8.5 Prepare SAM spatial query ------------
DefaultAssay(sobj_sam_ann) <- "SCT"
DefaultAssay(sobj_scrna_sam_ds) <- "SCT"

if (!"celltypes_v1" %in% colnames(sobj_sam_ann@meta.data)) sobj_sam_ann$celltypes_v1 <- sobj_sam_ann$celltypes

if ("predictions" %in% names(sobj_sam_ann@assays)) sobj_sam_ann[["predictions"]] <- NULL

old_pred_cols <- grep("^(predicted.id|prediction.score|scRNA_transfer_confidence|scRNA_predicted_celltype)", colnames(sobj_sam_ann@meta.data), value = TRUE)
if (length(old_pred_cols) > 0) sobj_sam_ann@meta.data[, old_pred_cols] <- NULL


## 8.15 Optional mixture-aware deconvolution with RCTD/spacexr ------------
## Seurat label transfer is not true deconvolution; RCTD estimates mixtures and is useful when one scRNA label maps to unexpected tissue zones.
library(spacexr)
library(SpatialExperiment)
library(SummarizedExperiment)

## Prepare matched raw-count matrices for spatial spots and the scRNA reference.
sam_spatial_counts <- GetAssayData(sobj_sam_ann, assay = "Spatial", layer = "counts")
sam_scrna_count_assay <- if ("RNA" %in% names(sobj_scrna_sam_ds@assays)) "RNA" else DefaultAssay(sobj_scrna_sam_ds)
sam_scrna_counts <- GetAssayData(sobj_scrna_sam_ds, assay = sam_scrna_count_assay, layer = "counts")
sam_shared_genes <- intersect(rownames(sam_spatial_counts), rownames(sam_scrna_counts))
sam_spatial_counts <- as.matrix(sam_spatial_counts[sam_shared_genes, ])
sam_scrna_counts <- as.matrix(sam_scrna_counts[sam_shared_genes, ])

## Build the spatial object from Seurat tissue coordinates for the current spacexr API.
sam_rctd_coords <- purrr::map_dfr(names(sobj_sam_ann@images), function(img) {
  co <- as.data.frame(GetTissueCoordinates(sobj_sam_ann, image = img, scale = NULL))
  x_col <- intersect(c("x", "imagecol", "col"), colnames(co))[1]; y_col <- intersect(c("y", "imagerow", "row"), colnames(co))[1]
  tibble(spot = rownames(co), x = co[[x_col]], y = co[[y_col]])
}) %>% distinct(spot, .keep_all = TRUE) %>% filter(spot %in% colnames(sam_spatial_counts))
sam_spatial_counts <- sam_spatial_counts[, sam_rctd_coords$spot, drop = FALSE]
sam_rctd_coords <- as.data.frame(sam_rctd_coords); rownames(sam_rctd_coords) <- sam_rctd_coords$spot
sam_rctd_coords <- sam_rctd_coords[colnames(sam_spatial_counts), , drop = FALSE]
sam_spatial_spe <- SpatialExperiment(assays = list(counts = sam_spatial_counts), spatialCoords = as.matrix(sam_rctd_coords[, c("x", "y")]))

## Build the scRNA reference and run RCTD in multi mode for mixed Visium-like spots.
sam_rctd_ref_coldata <- data.frame(cell_type = factor(sobj_scrna_sam_ds[[scrna_celltype_col]][colnames(sam_scrna_counts), 1]),
                                   row.names = colnames(sam_scrna_counts))
sam_reference_se <- SummarizedExperiment(assays = list(counts = sam_scrna_counts), colData = sam_rctd_ref_coldata)
sam_rctd_data <- createRctd(sam_spatial_spe, sam_reference_se, cell_type_col = "cell_type")
sam_rctd <- runRctd(sam_rctd_data, rctd_mode = "multi", max_cores = 1, max_multi_types = 4)
saveRDS(sam_rctd, file.path(sam_rctd_dir, "SAM_RCTD_multi_spacexr.rds"))

## Export full weights and the top three inferred cell types per spot.
sam_rctd_weight_assay <- if ("weights_full" %in% assayNames(sam_rctd)) "weights_full" else "weights"
sam_rctd_weights <- as.data.frame(t(as.matrix(assay(sam_rctd, sam_rctd_weight_assay)))) %>% rownames_to_column("spot")
write.csv(sam_rctd_weights, file.path(sam_rctd_table_dir, "SAM_RCTD_multi_celltype_weights.csv"), row.names = FALSE)
sam_rctd_top <- sam_rctd_weights %>% pivot_longer(-spot, names_to = "celltype", values_to = "weight") %>%
  group_by(spot) %>% arrange(desc(weight), .by_group = TRUE) %>% dplyr::slice_head(n = 3) %>% ungroup()
write.csv(sam_rctd_top, file.path(sam_rctd_table_dir, "SAM_RCTD_multi_top3_celltypes_per_spot.csv"), row.names = FALSE)

## Summarize global cell-type contribution to decide which RCTD labels are worth plotting.
sam_rctd_long <- sam_rctd_weights %>% pivot_longer(-spot, names_to = "celltype", values_to = "weight")
sam_rctd_celltype_summary <- sam_rctd_long %>% group_by(celltype) %>%
  summarise(mean_weight = mean(weight), max_weight = max(weight), n_spots_gt_0.1 = sum(weight > 0.1), n_spots_gt_0.2 = sum(weight > 0.2), .groups = "drop") %>%
  arrange(desc(mean_weight), desc(n_spots_gt_0.2))
write.csv(sam_rctd_celltype_summary, file.path(sam_rctd_table_dir, "SAM_RCTD_celltype_weight_summary.csv"), row.names = FALSE)

## Compare RCTD weights against curated spatial domains, similar to the Seurat prediction heatmap.
sam_rctd_ann_col <- if ("celltypes_v1" %in% colnames(sobj_sam_ann@meta.data)) "celltypes_v1" else "celltypes"
sam_rctd_ann <- sobj_sam_ann@meta.data %>% rownames_to_column("spot") %>% dplyr::select(spot, spatial_annotation = all_of(sam_rctd_ann_col))
sam_rctd_heat_df <- sam_rctd_long %>% left_join(sam_rctd_ann, by = "spot") %>% filter(!is.na(spatial_annotation)) %>%
  group_by(spatial_annotation, celltype) %>% summarise(mean_weight = mean(weight), .groups = "drop")
sam_rctd_heat_mat <- reshape2::acast(sam_rctd_heat_df, spatial_annotation ~ celltype, value.var = "mean_weight", fill = 0)
write.csv(sam_rctd_heat_mat, file.path(sam_rctd_table_dir, "SAM_RCTD_mean_weight_by_spatial_annotation.csv"))
sam_rctd_heat_mat_plot <- t(sam_rctd_heat_mat)
sam_rctd_heat_mat_plot <- sam_rctd_heat_mat_plot[order(rownames(sam_rctd_heat_mat_plot)), order(colnames(sam_rctd_heat_mat_plot)), drop = FALSE]
sam_rctd_heat_cap <- quantile(sam_rctd_heat_mat_plot[sam_rctd_heat_mat_plot > 0], 0.9, na.rm = TRUE)
sam_rctd_heat_cap <- max(0.05, as.numeric(sam_rctd_heat_cap))
col_fun_sam_rctd <- circlize::colorRamp2(c(0, sam_rctd_heat_cap / 2, sam_rctd_heat_cap), c("white", "orange", "maroon"))
ht_sam_rctd <- ComplexHeatmap::Heatmap(
  sam_rctd_heat_mat_plot, name = "Mean\nRCTD\nweight", col = col_fun_sam_rctd,
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  row_names_side = "left", column_names_side = "bottom", column_names_rot = 45,
  row_title = "scRNA-seq cell type", column_title = paste0("SAM spatial annotation: ", sam_rctd_ann_col),
  row_names_gp = grid::gpar(fontsize = 11), column_names_gp = grid::gpar(fontsize = 13),
  heatmap_legend_param = list(direction = "horizontal", title = paste0("Mean weight\n0-", signif(sam_rctd_heat_cap, 2), " cap"))
)
pdf(file.path(sam_rctd_plot_dir, "SAM_RCTD_mean_weight_by_spatial_annotation_heatmap.pdf"), width = 7, height = 8)
ComplexHeatmap::draw(ht_sam_rctd, heatmap_legend_side = "bottom")
dev.off()

## Add all RCTD weights to Seurat metadata so every scRNA cell type can be inspected spatially.
sam_rctd_celltypes <- setdiff(colnames(sam_rctd_weights), "spot")
sam_rctd_meta_features <- paste0("RCTD_", make.names(sam_rctd_celltypes))
for (i in seq_along(sam_rctd_celltypes)) sobj_sam_ann@meta.data[sam_rctd_weights$spot, sam_rctd_meta_features[i]] <- sam_rctd_weights[[sam_rctd_celltypes[i]]]

## Plot every RCTD cell-type weight in spatial maps; use pages so the file stays editable.
sam_rctd_spatial_images <- names(sobj_sam_ann@images)[intersect(c(1, 3), seq_along(names(sobj_sam_ann@images)))]
sam_rctd_feature_pages <- split(sam_rctd_meta_features, ceiling(seq_along(sam_rctd_meta_features) / 6))
pdf(file.path(sam_rctd_plot_dir, "SAM_RCTD_all_celltype_weights_spatial_images1_3.pdf"), width = 13, height = 16)
for (pg in seq_along(sam_rctd_feature_pages)) {
  print(SpatialFeaturePlot(sobj_sam_ann, features = sam_rctd_feature_pages[[pg]], images = sam_rctd_spatial_images,
                           crop = FALSE, ncol = 3, pt.size.factor = 1.5, alpha = c(0.1, 1)) +
          plot_annotation(title = paste0("SAM RCTD cell-type weights, page ", pg)))
}
dev.off()

## Plot dominant RCTD label and confidence to evaluate whether one label dominates each spot.
sam_rctd_dominant <- sam_rctd_long %>% group_by(spot) %>% dplyr::slice_max(weight, n = 1, with_ties = FALSE) %>%
  ungroup() %>% dplyr::rename(RCTD_dominant_celltype = celltype, RCTD_dominant_weight = weight)
sobj_sam_ann$RCTD_dominant_celltype <- NA_character_; sobj_sam_ann$RCTD_dominant_weight <- NA_real_
sobj_sam_ann@meta.data[sam_rctd_dominant$spot, "RCTD_dominant_celltype"] <- sam_rctd_dominant$RCTD_dominant_celltype
sobj_sam_ann@meta.data[sam_rctd_dominant$spot, "RCTD_dominant_weight"] <- sam_rctd_dominant$RCTD_dominant_weight
ggsave(file.path(sam_rctd_plot_dir, "SAM_RCTD_dominant_celltype_spatial_images1_3.pdf"),
       SpatialDimPlot(sobj_sam_ann, group.by = "RCTD_dominant_celltype", images = sam_rctd_spatial_images, crop = FALSE, ncol = 2, pt.size.factor = 1.6),
       width = 14, height = 7)
ggsave(file.path(sam_rctd_plot_dir, "SAM_RCTD_dominant_weight_spatial_images1_3.pdf"),
       SpatialFeaturePlot(sobj_sam_ann, features = "RCTD_dominant_weight", images = sam_rctd_spatial_images, crop = FALSE, ncol = 2, pt.size.factor = 1.6),
       width = 12, height = 6)

## Spatial pie map: show major RCTD mixture components per spot, plus Other for low-contribution cell types.
library(scatterpie) # install.packages("scatterpie") if this is missing
sam_rctd_coords_plot <- purrr::map_dfr(sam_rctd_spatial_images, function(img) {
  co <- as.data.frame(GetTissueCoordinates(sobj_sam_ann, image = img, scale = NULL))
  x_col <- intersect(c("x", "imagecol", "col"), colnames(co))[1]; y_col <- intersect(c("y", "imagerow", "row"), colnames(co))[1]
  tibble(spot = rownames(co), image = img, x = co[[x_col]], y = co[[y_col]])
}) %>% filter(spot %in% sam_rctd_weights$spot) %>%
  group_by(image) %>% mutate(x_plot = x - min(x, na.rm = TRUE), y_plot = y - min(y, na.rm = TRUE)) %>% ungroup()
sam_rctd_pie_celltypes <- sam_rctd_celltype_summary %>% filter(n_spots_gt_0.2 > 20) %>% dplyr::slice_head(n = 8) %>% pull(celltype)
if (length(sam_rctd_pie_celltypes) < 2) sam_rctd_pie_celltypes <- sam_rctd_celltype_summary %>% dplyr::slice_head(n = 8) %>% pull(celltype)
sam_rctd_pie_cols <- make.names(sam_rctd_pie_celltypes)
sam_rctd_pie_df <- sam_rctd_weights %>% dplyr::select(spot, all_of(sam_rctd_pie_celltypes)) %>% rename_with(make.names, -spot) %>%
  mutate(Other = pmax(0, 1 - rowSums(across(all_of(sam_rctd_pie_cols)), na.rm = TRUE)), pie_radius = 18) %>% right_join(sam_rctd_coords_plot, by = "spot")
sam_rctd_pie_cols <- c(sam_rctd_pie_cols, "Other")
pdf(file.path(sam_rctd_plot_dir, "SAM_RCTD_spatial_mixture_pie_images1_3.pdf"), width = 10, height = 8)
for (img in sam_rctd_spatial_images) {
  p_sam_rctd_pie <- ggplot(filter(sam_rctd_pie_df, image == img)) +
    geom_scatterpie(aes(x = x_plot, y = y_plot, r = pie_radius), cols = sam_rctd_pie_cols, color = NA, alpha = 0.9) +
    scale_fill_manual(values = c(scales::hue_pal()(length(sam_rctd_pie_celltypes)), Other = "grey80"),
                      labels = c(sam_rctd_pie_celltypes, "Other")) +
    scale_y_reverse() + coord_equal(expand = FALSE) + theme_classic() +
    labs(title = paste0("SAM RCTD spatial mixture pie map: ", img), fill = "scRNA cell type", x = NULL, y = NULL)
  print(p_sam_rctd_pie)
}
dev.off()

## Official spacexr visualizations: plotAllWeights is the RCTD/STdeconvolve-style pie chart.
library(cowplot) # used only to save the official plotAllWeights legend separately
sam_rctd_legend_spots <- sam_rctd_coords_plot$spot[sam_rctd_coords_plot$image == sam_rctd_spatial_images[1]]
sam_rctd_plotallweights_legend <- cowplot::get_legend(plotAllWeights(sam_rctd[, sam_rctd_legend_spots], assay_name = sam_rctd_weight_assay, r = 12, lwd = 0,
                                                                      title = "SAM RCTD mixtures") + theme(legend.position = "right"))
pdf(file.path(sam_rctd_plot_dir, "SAM_RCTD_plotAllWeights_official_color_legend.pdf"), width = 8, height = 10)
grid::grid.draw(sam_rctd_plotallweights_legend)
dev.off()

for (img in sam_rctd_spatial_images) {
  spots_i <- sam_rctd_coords_plot$spot[sam_rctd_coords_plot$image == img]
  img_width <- diff(range(sam_rctd_coords_plot$x[sam_rctd_coords_plot$image == img], na.rm = TRUE))
  img_height <- diff(range(sam_rctd_coords_plot$y[sam_rctd_coords_plot$image == img], na.rm = TRUE))
  pdf_w <- if (img_width >= img_height) 13 else 8; pdf_h <- if (img_width >= img_height) 5 else 10
  pdf(file.path(sam_rctd_plot_dir, paste0("SAM_RCTD_plotAllWeights_official_", img, ".pdf")), width = pdf_w, height = pdf_h)
  print(plotAllWeights(sam_rctd[, spots_i], assay_name = sam_rctd_weight_assay, r = 20, lwd = 0,
                       title = paste0("SAM RCTD mixtures: ", img)) + theme(legend.position = "none"))
  dev.off()
}

## Official spacexr single-cell-type weight maps for the strongest contributors.
sam_rctd_top_plot_celltypes <- sam_rctd_celltype_summary %>% dplyr::slice_head(n = 12) %>% pull(celltype)
sam_rctd_rotate_ccw_images <- c("sam_A_s1")
sam_rctd_flip_y_images <- c("sam_B_s1")
for (img in sam_rctd_spatial_images) {
  spots_i <- sam_rctd_coords_plot$spot[sam_rctd_coords_plot$image == img]
  sam_rctd_img <- sam_rctd[, spots_i]
  coords_i <- spatialCoords(sam_rctd_img)
  if (img %in% sam_rctd_rotate_ccw_images) coords_i <- cbind(x = -coords_i[, "y"], y = coords_i[, "x"])
  if (img %in% sam_rctd_flip_y_images) coords_i[, "y"] <- max(coords_i[, "y"], na.rm = TRUE) - coords_i[, "y"]
  spatialCoords(sam_rctd_img) <- coords_i
  pdf(file.path(sam_rctd_plot_dir, paste0("SAM_RCTD_plotCellTypeWeight_top_celltypes_official_", img, ".pdf")), width = 7, height = 7)
  for (ct in sam_rctd_top_plot_celltypes) print(plotCellTypeWeight(sam_rctd_img, cell_type = ct, assay_name = sam_rctd_weight_assay, size = 2, stroke = 0, title = paste0(ct, ": ", img)))
  dev.off()
}


  },
  trajectory = function() {
sobj_sam_ann <- qs::qread(atlas_input("saved_obj/sobj_sam_res0.4_final_v1_2026.7.24.qs"))
Seurat::Idents(sobj_sam_ann) <- sobj_sam_ann$celltypes
# 9 Monocle3 trajectory analysis for SAM ####
# --------------------------------------- #

library(monocle3)
library(SeuratWrappers)
library(qs)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(clusterProfiler)
library(forcats)
library(stringr)

dir.create("trajectory/SAM/monocle3", showWarnings = FALSE, recursive = TRUE)
dir.create("tables/trajectory", showWarnings = FALSE, recursive = TRUE)
sam_pseudotime_colors <- viridisLite::plasma(100)
sam_pseudotime_anno_colors <- list(Pseudotime = sam_pseudotime_colors)

## 9.1 Prepare SAM trajectory input ------------

# Input = re-cleaned SAM object (6 coarse domains in celltypes_v1) with NMF k=16 factor scores
# and the trichome score already attached as per-spot meta.data columns.
if (!exists("sobj_sam_ann")) sobj_sam_ann <- qread("saved_obj/sobj_sam_split_cleaned_res0.4_traits_2026.7.9.qs")

DefaultAssay(sobj_sam_ann) <- "SCT"
Idents(sobj_sam_ann) <- sobj_sam_ann$celltypes_v1

# Full SAM trajectory includes all six newly annotated SAM domains.
sam <- sobj_sam_ann
sam$sam_traj_celltype <- as.character(sam$celltypes_v1)

table(sam$sam_traj_celltype, useNA = "ifany")

pdf("trajectory/SAM/monocle3/SAM_monocle3_input_check.pdf", width = 10, height = 9)
print(
  DimPlot(sam, reduction = "umap", group.by = "sam_traj_celltype",
          label = TRUE, repel = TRUE) +
    theme_classic()
)
print(SpatialDimPlot(sam, group.by = "sam_traj_celltype", crop = FALSE, ncol = 3))
dev.off()

## 9.2 Convert SAM Seurat object to Monocle3 cell_data_set ------------

cds_sam <- as.cell_data_set(sam)

fData(cds_sam)$gene_short_name <- rownames(fData(cds_sam))
colData(cds_sam)$sam_traj_celltype <- sam$sam_traj_celltype[match(colnames(cds_sam), colnames(sam))]
colData(cds_sam)$seurat_cluster <- as.character(sam$seurat_clusters[match(colnames(cds_sam), colnames(sam))])

reducedDims(cds_sam)$UMAP <- Embeddings(sam, reduction = "umap")[colnames(cds_sam), ]


## 9.3 Cluster cells and learn full SAM Monocle3 graph ------------

sam_monocle_resolution <- 0.002
sam_minimal_branch_len <- 12

cds_sam <- cluster_cells(cds_sam, reduction_method = "UMAP", resolution = sam_monocle_resolution)

colData(cds_sam)$monocle3_cluster <- as.character(clusters(cds_sam))
sam$monocle3_cluster <- colData(cds_sam)$monocle3_cluster[match(colnames(sam), colnames(cds_sam))]

p_sam_celltype <- plot_cells(
  cds_sam, color_cells_by = "sam_traj_celltype", show_trajectory_graph = FALSE,
  label_cell_groups = TRUE, label_groups_by_cluster = FALSE,
  cell_size = 0.8, group_label_size = 4
) + ggtitle("SAM Monocle3 input colored by annotation")

ggsave("trajectory/SAM/monocle3/SAM_monocle3_input_celltypes.pdf", p_sam_celltype, width = 6, height = 6)

p_sam_monocle_cluster <- plot_cells(
  cds_sam, color_cells_by = "monocle3_cluster", show_trajectory_graph = FALSE,
  label_cell_groups = TRUE, label_groups_by_cluster = TRUE,
  cell_size = 0.8, group_label_size = 4
) + ggtitle(paste0("SAM Monocle3 clusters; resolution = ", sam_monocle_resolution))

ggsave("trajectory/SAM/monocle3/SAM_monocle3_clusters.pdf", p_sam_monocle_cluster, width = 6, height = 6)

sam_monocle_celltype_table <- table(colData(cds_sam)$monocle3_cluster, colData(cds_sam)$sam_traj_celltype)
sam_monocle_seurat_table <- table(colData(cds_sam)$monocle3_cluster, colData(cds_sam)$seurat_cluster)

write.csv(as.data.frame.matrix(sam_monocle_celltype_table),
          "tables/trajectory/SAM_monocle3_cluster_vs_celltype.csv")

write.csv(as.data.frame.matrix(sam_monocle_seurat_table),
          "tables/trajectory/SAM_monocle3_cluster_vs_seurat_cluster.csv")

cds_sam <- learn_graph(
  cds_sam, use_partition = FALSE, close_loop = FALSE,
  learn_graph_control = list(prune_graph = TRUE, minimal_branch_len = sam_minimal_branch_len)
)

p_sam_graph_raw <- plot_cells(
  cds_sam, color_cells_by = "sam_traj_celltype",
  label_cell_groups = TRUE, label_groups_by_cluster = TRUE,
  label_branch_points = FALSE, label_roots = FALSE,
  label_leaves = FALSE, cell_size = 1, group_label_size = 4
) + ggtitle("SAM Monocle3 graph before cluster recoding")

ggsave("trajectory/SAM/monocle3/SAM_monocle3_full_graph_raw.pdf", p_sam_graph_raw, width = 8, height = 7)

## 9.4 Highlight individual Monocle3 clusters ------------

sam_clst_ids <- sort(unique(as.character(clusters(cds_sam))))

pdf("trajectory/SAM/monocle3/SAM_monocle3_subclusters_check_UMAP.pdf", width = 7, height = 7)
for (i in sam_clst_ids) {
  cells_i <- names(clusters(cds_sam)[as.character(clusters(cds_sam)) == i])
  print(
    DimPlot(sam, reduction = "umap", cells.highlight = cells_i,
            pt.size = 0.5, sizes.highlight = 1.2) +
      ggtitle(paste0("SAM Monocle3 cluster ", i, " | n = ", length(cells_i)))
      )
}
dev.off()

pdf("trajectory/SAM/monocle3/SAM_monocle3_subclusters_check_spatial.pdf", width = 13, height = 8)
for (i in sam_clst_ids) {
  cells_i <- names(clusters(cds_sam)[as.character(clusters(cds_sam)) == i])
  print(
    SpatialDimPlot(
      sam, cells.highlight = cells_i, alpha = c(0.6, 1),
      pt.size.factor = 7, crop = TRUE, images = c('sam_B_s1','sam_B_s2')
    ) + ggtitle(paste0("SAM Monocle3 cluster ", i, " | n = ", length(cells_i))))
  }
dev.off()


## 9.5 Define SAM trajectory labels ------------

# Use the new Seurat annotation directly; Monocle3 clusters remain available
# separately for fine-resolution topology and trajectory fitting.
colData(cds_sam)$assigned_cell_type <- as.character(colData(cds_sam)$sam_traj_celltype)

sam$assigned_cell_type_monocle3 <- colData(cds_sam)$assigned_cell_type[match(colnames(sam), colnames(cds_sam))]

sobj_sam_ann$assigned_cell_type_monocle3 <- NA
sobj_sam_ann@meta.data[colnames(sam), "assigned_cell_type_monocle3"] <- sam$assigned_cell_type_monocle3

assigned_celltype_table <- table(colData(cds_sam)$assigned_cell_type, colData(cds_sam)$sam_traj_celltype)
assigned_seurat_table <- table(colData(cds_sam)$assigned_cell_type, colData(cds_sam)$seurat_cluster)

write.csv(as.data.frame.matrix(assigned_celltype_table),
          "tables/trajectory/SAM_monocle3_assigned_celltype_vs_original_celltype.csv")

write.csv(as.data.frame.matrix(assigned_seurat_table),
          "tables/trajectory/SAM_monocle3_assigned_celltype_vs_seurat_cluster.csv")

p_sam_assigned <- plot_cells(
  cds_sam, color_cells_by = "assigned_cell_type",
  label_cell_groups = TRUE, label_groups_by_cluster = FALSE,
  label_branch_points = FALSE, label_roots = FALSE,
  label_leaves = FALSE, cell_size = 1, group_label_size = 4
) + ggtitle("SAM Monocle3 assigned cell types")

ggsave("trajectory/SAM/monocle3/SAM_monocle3_assigned_celltypes.pdf", p_sam_assigned, width = 8, height = 7)

## 9.6 Order full SAM graph in pseudotime as a topology diagnostic ------------

## Full SAM contains multiple developmental branches, so pseudotime distances are not directly comparable across all leaves, pith, and vascular cells.
## Use branch-specific pseudotime below for final ordered heatmaps and gene-trend interpretation.

root_label_sam <- "Meristem"

if (!root_label_sam %in% colData(cds_sam)$assigned_cell_type) {
  stop("Meristem was not found in assigned_cell_type. Check SAM annotation labels.")
}

root_cells_sam <- colnames(cds_sam)[colData(cds_sam)$assigned_cell_type == root_label_sam]

cds_sam <- order_cells(cds_sam, reduction_method = "UMAP", root_cells = root_cells_sam)

cds_sam$monocle3_pseudotime <- pseudotime(cds_sam)
sam$monocle3_pseudotime <- pseudotime(cds_sam)[match(colnames(sam), colnames(cds_sam))]

sam_pt_range <- range(sam$monocle3_pseudotime, na.rm = TRUE)

p_sam_pt <- plot_cells(
  cds_sam, color_cells_by = "pseudotime",
  label_cell_groups = TRUE, label_groups_by_cluster = FALSE,
  label_branch_points = FALSE, label_roots = FALSE,
  label_leaves = FALSE, cell_size = 1.1
) +
  scale_color_gradientn(colors = sam_pseudotime_colors,
                        limits = sam_pt_range, name = "Pseudotime") +
  ggtitle("SAM full-tissue Monocle3 pseudotime diagnostic; root = Meristem")

ggsave("trajectory/SAM/monocle3/SAM_monocle3_full_pseudotime.pdf",p_sam_pt, width = 8, height = 7)


## 9.7 Export full SAM pseudotime diagnostic to Seurat ------------

sobj_sam_ann$monocle3_pseudotime <- NA_real_
sobj_sam_ann@meta.data[colnames(sam), "monocle3_pseudotime"] <- sam$monocle3_pseudotime

sam_pt_meta <- as.data.frame(colData(cds_sam)) %>%
  rownames_to_column("spot") %>%
  filter(!is.na(monocle3_pseudotime))

write.csv(sam_pt_meta, "tables/trajectory/SAM_monocle3_full_pseudotime_metadata.csv", row.names = FALSE)

sam_celltype_order <- sam_pt_meta %>%
  group_by(assigned_cell_type) %>%
  summarize(median_pseudotime = median(monocle3_pseudotime, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_pseudotime) %>%
  pull(assigned_cell_type)

sam_pt_meta$assigned_cell_type <- factor(sam_pt_meta$assigned_cell_type,
                                         levels = rev(sam_celltype_order))

p_sam_pt_box <- ggplot(sam_pt_meta, aes(x = monocle3_pseudotime, y = assigned_cell_type,
                                        fill = assigned_cell_type)) +
  geom_boxplot(outlier.size = 0.25) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Monocle3 pseudotime", y = NULL)

ggsave("trajectory/SAM/monocle3/SAM_monocle3_full_pseudotime_by_celltype.pdf", p_sam_pt_box, width = 8, height = 6)

p_sam_pt_spatial <- SpatialFeaturePlot(sam, features = "monocle3_pseudotime", crop = FALSE, ncol = 3) &
  scale_fill_gradientn(colors = sam_pseudotime_colors, limits = sam_pt_range, name = "Pseudotime", na.value = "grey90") &
  theme(legend.position = "right")
ggsave("trajectory/SAM/monocle3/SAM_monocle3_full_pseudotime_spatial.pdf", p_sam_pt_spatial, width = 14, height = 9)


## 9.12 Upward SAM branch: Meristem to leaf primordium ------------

# Add "Epidermis" below only if the initial branch graph supports its inclusion.
sam_upward_types <- c("Meristem", "Leaf primordium")

sam_upward_types <- intersect(sam_upward_types, unique(colData(cds_sam)$assigned_cell_type))

cds_sam_up <- cds_sam[, colData(cds_sam)$assigned_cell_type %in% sam_upward_types]

sam_up_resolution <- 0.03
sam_up_minimal_branch_len <- 10

cds_sam_up <- cluster_cells(cds_sam_up, reduction_method = "UMAP", resolution = sam_up_resolution)

cds_sam_up <- learn_graph(
  cds_sam_up, use_partition = FALSE, close_loop = FALSE,
  learn_graph_control = list(prune_graph = TRUE, minimal_branch_len = sam_up_minimal_branch_len))

root_cells_sam_up <- colnames(cds_sam_up)[colData(cds_sam_up)$assigned_cell_type == "Meristem"]

cds_sam_up <- order_cells(cds_sam_up, reduction_method = "UMAP", root_cells = root_cells_sam_up)
cds_sam_up$branch_pseudotime <- pseudotime(cds_sam_up)

sam$monocle3_pseudotime_upward_leaf <- NA_real_
sam@meta.data[colnames(cds_sam_up), "monocle3_pseudotime_upward_leaf"] <- cds_sam_up$branch_pseudotime
sam$monocle3_cluster_upward_leaf <- NA_character_
sam@meta.data[colnames(cds_sam_up), "monocle3_cluster_upward_leaf"] <- as.character(clusters(cds_sam_up))

sobj_sam_ann$monocle3_pseudotime_upward_leaf <- NA_real_
sobj_sam_ann@meta.data[colnames(cds_sam_up), "monocle3_pseudotime_upward_leaf"] <- cds_sam_up$branch_pseudotime
sobj_sam_ann$monocle3_cluster_upward_leaf <- NA_character_
sobj_sam_ann@meta.data[colnames(cds_sam_up), "monocle3_cluster_upward_leaf"] <- as.character(clusters(cds_sam_up))

## 9.12.1 Rename upward-branch Monocle3 clusters after checking histology ------------
## Start from Monocle3 cluster IDs, then manually assign biological subcluster labels.
colData(cds_sam_up)$monocle3_cluster_upward_leaf <- as.character(clusters(cds_sam_up))
colData(cds_sam_up)$monocle3_subcluster_upward_leaf <- as.character(colData(cds_sam_up)$monocle3_cluster_upward_leaf)

colData(cds_sam_up)$monocle3_subcluster_upward_leaf <- dplyr::recode(
  colData(cds_sam_up)$monocle3_subcluster_upward_leaf,
  "7" = "Meristem",
  '15' = "Leaf primordium2",'16' = "Leaf primordium2",
  "8" = "Leaf primordium1",'12' = "Leaf primordium1","9" = "Leaf primordium1",'13' = "Leaf primordium1",
  "3" = "Leaf primordium3 ","4" = "Leaf primordium3","5" = "Leaf primordium3",
  "1" = "Leaf primordium4","6" = "Leaf primordium4",
  "10" = "Leaf primordium5",'11' = "Leaf primordium5",
  '14' = "Leaf primordium6",'17' = "Leaf primordium6","2" = "Leaf primordium6",
  .default = colData(cds_sam_up)$monocle3_subcluster_upward_leaf
)
colData(cds_sam_up)$monocle3_subcluster_upward_leaf <- gsub("\\s+", " ", trimws(as.character(colData(cds_sam_up)$monocle3_subcluster_upward_leaf)))

dir.create("tables/trajectory", showWarnings = FALSE, recursive = TRUE)
sam_up_cluster_label_map <- as.data.frame(colData(cds_sam_up)) %>% rownames_to_column("spot") %>%
  group_by(monocle3_cluster_upward_leaf, monocle3_subcluster_upward_leaf) %>%
  summarize(n_spots = n(), dominant_original_celltype = names(sort(table(assigned_cell_type), decreasing = TRUE))[1], .groups = "drop") %>%
  arrange(as.numeric(monocle3_cluster_upward_leaf))
write.csv(sam_up_cluster_label_map, "tables/trajectory/SAM_upward_leaf_branch_monocle3_subcluster_labels.csv", row.names = FALSE)
sam$monocle3_subcluster_upward_leaf <- NA_character_
sam@meta.data[colnames(cds_sam_up), "monocle3_subcluster_upward_leaf"] <- colData(cds_sam_up)$monocle3_subcluster_upward_leaf
sobj_sam_ann$monocle3_subcluster_upward_leaf <- NA_character_
sobj_sam_ann@meta.data[colnames(cds_sam_up), "monocle3_subcluster_upward_leaf"] <- colData(cds_sam_up)$monocle3_subcluster_upward_leaf
sam_up_monocle3_meta <- tibble(spot = colnames(cds_sam_up),
                               monocle3_pseudotime_upward_leaf = as.numeric(cds_sam_up$branch_pseudotime),
                               monocle3_cluster_upward_leaf = as.character(colData(cds_sam_up)$monocle3_cluster_upward_leaf),
                               monocle3_subcluster_upward_leaf = as.character(colData(cds_sam_up)$monocle3_subcluster_upward_leaf))
write.csv(sam_up_monocle3_meta, "tables/trajectory/SAM_upward_leaf_branch_monocle3_metadata.csv", row.names = FALSE)

## First-panel QC for the SAM upward figure: original Seurat labels vs renamed Monocle3 subclusters.
sam_up_figure_images <- names(sobj_sam_ann@images)[intersect(c(1, 3,4), seq_along(names(sobj_sam_ann@images)))]
sam_up_spatial_obj <- subset(sobj_sam_ann, cells = colnames(cds_sam_up))
sam_up_original_label_col <- if ("celltypes_v1" %in% colnames(sam_up_spatial_obj@meta.data)) "celltypes_v1" else "celltypes"
p_sam_up_original_spatial <- SpatialDimPlot(sam_up_spatial_obj, group.by = sam_up_original_label_col, images = sam_up_figure_images, crop = FALSE, ncol = 2, pt.size.factor = 1.8) +
  plot_annotation(title = "SAM upward branch: original Seurat labels")
p_sam_up_monocle_spatial <- SpatialDimPlot(sam_up_spatial_obj, group.by = "monocle3_subcluster_upward_leaf", images = sam_up_figure_images, crop = FALSE, ncol = 2, pt.size.factor = 1.8) +
  plot_annotation(title = "SAM upward branch: renamed Monocle3 subclusters")
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_original_seurat_labels_spatial_images1_3_4.pdf", p_sam_up_original_spatial, width = 10, height = 5)
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_renamed_monocle3_subclusters_spatial_images1_3_4.pdf", p_sam_up_monocle_spatial, width = 10, height = 5)

# qsave(sam, "saved_obj/sobj_sam_split_cleaned_res0.4_monocle3_upward_subclusters.qs")
qsave(sobj_sam_ann, "saved_obj/sobj_sam_ann_monocle3_upward_subclusters.2026.7.15.qs")

## Use the same purple-to-yellow pseudotime scale if this section is rerun alone.
if (!exists("sam_pseudotime_colors")) sam_pseudotime_colors <- viridisLite::plasma(100)
if (!exists("sam_pseudotime_anno_colors")) sam_pseudotime_anno_colors <- list(Pseudotime = sam_pseudotime_colors)

p_sam_up_monocle <- plot_cells(
  cds_sam_up, color_cells_by = "monocle3_subcluster_upward_leaf",
  label_cell_groups = TRUE, label_groups_by_cluster = FALSE,
  label_branch_points = FALSE, label_roots = FALSE,
  label_leaves = FALSE, cell_size = 1.1
) +
  ggtitle("SAM upward branch Monocle3 UMAP: renamed subclusters")

ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_renamed_monocle3_subclusters_UMAP.pdf", p_sam_up_monocle, width = 6, height = 5)

p_sam_up_pt <- plot_cells(
  cds_sam_up, color_cells_by = "pseudotime",
  label_cell_groups = FALSE, label_groups_by_cluster = FALSE,
  label_branch_points = FALSE, label_roots = FALSE,
  label_leaves = FALSE, cell_size = 1.1
) +
  scale_color_gradientn(colors = sam_pseudotime_colors, name = "Pseudotime") +
  ggtitle("SAM upward branch: Meristem to leaf primordium")

ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_pseudotime.pdf", p_sam_up_pt, width = 7, height = 5)

## Boxplots for upward pseudotime by original Seurat labels and renamed Monocle3 subclusters.
sam_up_pt_box_df <- as.data.frame(colData(cds_sam_up)) %>% rownames_to_column("spot") %>%
  mutate(branch_pseudotime = as.numeric(branch_pseudotime))
sam_up_seurat_pt_order <- sam_up_pt_box_df %>% group_by(assigned_cell_type) %>%
  summarize(pt_median = median(branch_pseudotime, na.rm = TRUE), .groups = "drop") %>% arrange(pt_median) %>% pull(assigned_cell_type)
sam_up_monocle_pt_order <- sam_up_pt_box_df %>% group_by(monocle3_subcluster_upward_leaf) %>%
  summarize(pt_median = median(branch_pseudotime, na.rm = TRUE), .groups = "drop") %>% arrange(pt_median) %>% pull(monocle3_subcluster_upward_leaf)
sam_up_pt_box_df <- sam_up_pt_box_df %>%
  mutate(assigned_cell_type = factor(assigned_cell_type, levels = sam_up_seurat_pt_order),
         monocle3_subcluster_upward_leaf = factor(monocle3_subcluster_upward_leaf, levels = sam_up_monocle_pt_order))
p_sam_up_pt_box_seurat <- ggplot(sam_up_pt_box_df, aes(assigned_cell_type, branch_pseudotime, fill = assigned_cell_type)) +
  geom_boxplot(outlier.size = 0.25) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "SAM upward pseudotime by original Seurat label", x = NULL, y = "Upward branch pseudotime")
p_sam_up_pt_box_monocle <- ggplot(sam_up_pt_box_df, aes(monocle3_subcluster_upward_leaf, branch_pseudotime, fill = monocle3_subcluster_upward_leaf)) +
  geom_boxplot(outlier.size = 0.25) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "SAM upward pseudotime by renamed Monocle3 subcluster", x = NULL, y = "Upward branch pseudotime")
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_pseudotime_boxplot_original_seurat_labels.pdf", p_sam_up_pt_box_seurat, width = 7, height = 5)
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_pseudotime_boxplot_renamed_monocle3_subclusters.pdf", p_sam_up_pt_box_monocle, width = 9, height = 5)

## Check the upward branch by Monocle3 subcluster and by spatial position.
p_sam_up_cluster <- plot_cells(
  cds_sam_up, color_cells_by = "cluster",
  label_cell_groups = TRUE, label_groups_by_cluster = TRUE,
  label_branch_points = FALSE, label_roots = FALSE,
  label_leaves = FALSE, cell_size = 1.1
) + ggtitle("SAM upward branch Monocle3 subclusters")
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_monocle3_clusters.pdf", p_sam_up_cluster, width = 7, height = 5)

p_sam_up_subcluster <- plot_cells(
  cds_sam_up, color_cells_by = "monocle3_subcluster_upward_leaf",
  label_cell_groups = TRUE, label_groups_by_cluster = FALSE,
  label_branch_points = FALSE, label_roots = FALSE,
  label_leaves = FALSE, cell_size = 1.1
) + ggtitle("SAM upward branch curated Monocle3 subclusters")
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_curated_monocle3_subclusters.pdf", p_sam_up_subcluster, width = 7, height = 5)

## Highlight each upward-branch Monocle3 cluster separately for manual topology inspection.
sam_up_clst_ids <- sort(unique(as.character(clusters(cds_sam_up))))
pdf("trajectory/SAM/monocle3/SAM_upward_leaf_branch_monocle3_subclusters_check_UMAP.pdf", width = 7, height = 7)
for (i in sam_up_clst_ids) {
  cells_i <- names(clusters(cds_sam_up)[as.character(clusters(cds_sam_up)) == i])
  print(DimPlot(sam, reduction = "umap", cells.highlight = cells_i, pt.size = 0.5, sizes.highlight = 1.2) +
          ggtitle(paste0("SAM upward branch Monocle3 cluster ", i, " | n = ", length(cells_i))))
}
dev.off()
sam_up_spatial_images <- names(sam@images)[c(1, 3)]
sam_up_spatial_images <- sam_up_spatial_images[!is.na(sam_up_spatial_images)]
pdf("trajectory/SAM/monocle3/SAM_upward_leaf_branch_monocle3_subclusters_check_spatial.pdf", width = 12, height = 6)
for (i in sam_up_clst_ids) {
  cells_i <- names(clusters(cds_sam_up)[as.character(clusters(cds_sam_up)) == i])
  print(SpatialDimPlot(sam, cells.highlight = cells_i, images = sam_up_spatial_images,
                       alpha = c(0.5, 1), pt.size.factor = 2, crop = FALSE, ncol = 2) +
          ggtitle(paste0("SAM upward branch Monocle3 cluster ", i, " | n = ", length(cells_i))))
}
dev.off()

## Color the same upward branch by original Seurat cell-type annotation copied from celltypes_v1.
## cell type names available in: colnames(colData(cds_sam_up))
p_sam_up_seurat_celltype <- plot_cells(
  cds_sam_up, color_cells_by = "assigned_cell_type",
  label_cell_groups = TRUE, label_groups_by_cluster = FALSE,
  label_branch_points = FALSE, label_roots = FALSE,
  label_leaves = FALSE, cell_size = 1.1
) + ggtitle("SAM upward branch original Seurat cell types")
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_original_seurat_celltypes.pdf", p_sam_up_seurat_celltype, width = 7, height = 5)

p_sam_up_pt_spatial <- SpatialFeaturePlot(sam, features = "monocle3_pseudotime_upward_leaf", crop = FALSE, ncol = 3) &
  scale_fill_gradientn(colors = sam_pseudotime_colors, name = "Pseudotime", na.value = "grey90") &
  theme(legend.position = "right")
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_pseudotime_spatial.pdf", p_sam_up_pt_spatial, width = 14, height = 9)

p_sam_up_seurat_celltype_spatial <- SpatialDimPlot(sam[, colnames(cds_sam_up)], group.by = "sam_traj_celltype", crop = FALSE, ncol = 3) +
  plot_annotation(title = "SAM upward branch original Seurat cell types in spatial sections")
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_original_seurat_celltypes_spatial.pdf", p_sam_up_seurat_celltype_spatial, width = 14, height = 9)

p_sam_up_cluster_spatial <- SpatialDimPlot(sam, group.by = "monocle3_cluster_upward_leaf", crop = FALSE, ncol = 3) +
  plot_annotation(title = "SAM upward branch Monocle3 subclusters in spatial sections")
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_monocle3_clusters_spatial.pdf", p_sam_up_cluster_spatial, width = 14, height = 14)

p_sam_up_curated_subcluster_spatial <- SpatialDimPlot(sam, group.by = "monocle3_subcluster_upward_leaf", crop = FALSE, ncol = 3) +
  plot_annotation(title = "SAM upward branch curated Monocle3 subclusters in spatial sections")
ggsave("trajectory/SAM/monocle3/SAM_upward_leaf_branch_curated_monocle3_subclusters_spatial.pdf", p_sam_up_curated_subcluster_spatial, width = 14, height = 11)

sam_up_cluster_levels <- sort(unique(na.omit(sam$monocle3_cluster_upward_leaf)))
pdf("trajectory/SAM/monocle3/SAM_upward_leaf_branch_monocle3_cluster_highlight_spatial.pdf", width = 14, height = 9)
for (cluster_i in sam_up_cluster_levels) {
    cluster_label_i <- paste0("Cluster_", cluster_i)
    sam$upward_cluster_highlight <- ifelse(sam$monocle3_cluster_upward_leaf == cluster_i, cluster_label_i, "Other")
    sam$upward_cluster_highlight[is.na(sam$upward_cluster_highlight)] <- "Not_in_upward_branch"
    cluster_cols_i <- c(setNames("#D95F02", cluster_label_i), "Not_in_upward_branch" = "grey90", "Other" = "grey70")
    print(SpatialDimPlot(sam, group.by = "upward_cluster_highlight", crop = FALSE, ncol = 3, cols = cluster_cols_i) +
            plot_annotation(title = paste0("SAM upward branch Monocle3 cluster ", cluster_i)))
}
dev.off()
sam$upward_cluster_highlight <- NULL

## Dynamic-gene cutoffs; keep local defaults so the upward block can be rerun alone.
if (!exists("sam_monocle_dynamic_q_cutoff")) sam_monocle_dynamic_q_cutoff <- 0.01
if (!exists("sam_monocle_dynamic_morans_cutoff")) sam_monocle_dynamic_morans_cutoff <- 0
if (!exists("sam_monocle_dynamic_filter_label")) sam_monocle_dynamic_filter_label <- paste0("Monocle3 graph_test: q_value < ", sam_monocle_dynamic_q_cutoff, ", Moran's I > ", sam_monocle_dynamic_morans_cutoff)

deg_sam_up <- graph_test(cds_sam_up, neighbor_graph = "principal_graph", cores = 4)

deg_sam_up$gene <- rownames(deg_sam_up)
deg_sam_up <- deg_sam_up %>% arrange(q_value, desc(morans_I))

sam_up_genes <- deg_sam_up %>%
  filter(q_value < sam_monocle_dynamic_q_cutoff, morans_I > sam_monocle_dynamic_morans_cutoff) %>%

qs::qsave(cds_sam_up, "trajectory/SAM/monocle3/SAM_upward_leaf_branch_monocle3_cds.qs")
qs::qsave(sobj_sam_ann, "saved_obj/sobj_sam_res0.4_final_v1_2026.7.24.qs")

  },
  wgcna = function() {
sobj_sam_ann <- qs::qread(atlas_input("saved_obj/sobj_sam_res0.4_final_v1_2026.7.24.qs"))
Seurat::Idents(sobj_sam_ann) <- sobj_sam_ann$celltypes
# 11 Final hdWGCNA after annotation, NMF, and trajectory analysis ####
# --------------------------------------- #

## Final hdWGCNA is run LAST, on the developmental core only (Meristem + Leaf primordium =
## 1,426 spots), so that the co-expression network is learned on the tissue where the trajectory
## and NMF programs are active, and so it can be correlated against the finished trajectory traits.
library(hdWGCNA); library(WGCNA); library(UCell); library(harmony)
library(dplyr); library(ggplot2); library(patchwork); library(tibble); library(tidyr)
library(ComplexHeatmap); library(circlize); library(clusterProfiler)
dir.create("hdWGCNA/SAM", showWarnings = FALSE, recursive = TRUE)
dir.create("hdWGCNA/SAM/TOM", showWarnings = FALSE, recursive = TRUE)
dir.create("tables/hdWGCNA", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/SAM", showWarnings = FALSE, recursive = TRUE)
allowWGCNAThreads(); set.seed(12345)

## 11.1 Subset to the developmental core ------------
# celltypes_v1 contains NA for a few spots, so subset by explicit cell names (subset(subset=...)
# on a column with NA drops the wrong cells). Meristem + Leaf primordium are the trajectory core.
if (!exists("sobj_sam_ann")) sobj_sam_ann <- qread("saved_obj/sobj_sam_split_cleaned_res0.4_monocle3_2026.7.9.qs")
DefaultAssay(sobj_sam_ann) <- "SCT"
sam_dev_types <- c("Meristem", "Leaf primordium")
md <- sobj_sam_ann@meta.data

## Build a stable Meristem/Leaf primordium label for hdWGCNA, using the best available annotation.
sam_hdwgcna_label_source <- if (all(sam_dev_types %in% unique(na.omit(md$celltypes_v1)))) "celltypes_v1" else if (all(sam_dev_types %in% unique(na.omit(md$assigned_cell_type_monocle3)))) "assigned_cell_type_monocle3" else if ("monocle3_subcluster_upward_leaf" %in% colnames(md)) "monocle3_subcluster_upward_leaf" else NA_character_
if (is.na(sam_hdwgcna_label_source)) stop("No annotation column contains Meristem and Leaf primordium for hdWGCNA.")
sobj_sam_ann$sam_dev_hdwgcna_group <- as.character(sobj_sam_ann@meta.data[[sam_hdwgcna_label_source]])
if (sam_hdwgcna_label_source == "monocle3_subcluster_upward_leaf") {
  sobj_sam_ann$sam_dev_hdwgcna_group <- case_when(
    sobj_sam_ann$sam_dev_hdwgcna_group == "Meristem" ~ "Meristem",
    grepl("Leaf primordium", sobj_sam_ann$sam_dev_hdwgcna_group) ~ "Leaf primordium",
    TRUE ~ NA_character_
  )
}
print(table(sobj_sam_ann$sam_dev_hdwgcna_group, useNA = "ifany"))
sam_dev_cells <- rownames(sobj_sam_ann@meta.data)[!is.na(sobj_sam_ann$sam_dev_hdwgcna_group) & sobj_sam_ann$sam_dev_hdwgcna_group %in% sam_dev_types]
if (length(sam_dev_cells) == 0) stop("No SAM developmental cells found for hdWGCNA.")
sub <- subset(sobj_sam_ann, cells = sam_dev_cells)
sub$sam_dev_hdwgcna_group <- droplevels(factor(sub$sam_dev_hdwgcna_group, levels = sam_dev_types))
sub$hdWGCNA_group <- "SAM_dev_core"
sam_dev_types <- intersect(sam_dev_types, levels(sub$sam_dev_hdwgcna_group))
if (length(sam_dev_types) < 2) stop("hdWGCNA needs both Meristem and Leaf primordium. Found: ", paste(levels(sub$sam_dev_hdwgcna_group), collapse = ", "))

## 11.2 Set up WGCNA gene selection ------------
# Keep genes expressed in >=5% of spots; store network under name "SAM_dev".
sub <- SetupForWGCNA(sub, gene_select = "fraction", fraction = 0.05, wgcna_name = "SAM_dev")

## 11.3 Construct metacells ------------
# Aggregate similar spots into metacells within each domain x section to denoise the network;
# k=15 nearest neighbors, max_shared=8 caps metacell overlap, min_cells=20 avoids dropping small Meristem groups.
sub <- MetacellsByGroups(sub, group.by = c("sam_dev_hdwgcna_group", "orig.ident", "hdWGCNA_group"), ident.group = "sam_dev_hdwgcna_group",
                         reduction = "pca", k = 15, max_shared = 8, min_cells = 20, wgcna_name = "SAM_dev")
sub <- NormalizeMetacells(sub)
print(table(sub$sam_dev_hdwgcna_group, useNA = "ifany"))
print(table(GetMetacellObject(sub, wgcna_name = "SAM_dev")$hdWGCNA_group, useNA = "ifany"))

## 11.4 Set expression matrix and choose soft-thresholding power ------------
sam_dev_types <- intersect(sam_dev_types, unique(as.character(sub$sam_dev_hdwgcna_group)))
if (length(sam_dev_types) < 2) stop("Metacell construction dropped one developmental group. Found: ", paste(unique(as.character(sub$sam_dev_hdwgcna_group)), collapse = ", "))
sub <- SetDatExpr(sub, group_name = "SAM_dev_core", group.by = "hdWGCNA_group", assay = "SCT", slot = "data", wgcna_name = "SAM_dev")
sub <- TestSoftPowers(sub, networkType = "signed")

pdf("hdWGCNA/SAM/SAM_dev_soft_power_diagnostics.pdf", width = 10, height = 8)
print(wrap_plots(PlotSoftPowers(sub), ncol = 2))
dev.off()

# Pick the first power reaching scale-free topology fit R^2 >= 0.8 (here soft power = 8).
sft <- GetPowerTable(sub)
sam_soft_power <- sft$Power[which(sft$SFT.R.sq >= 0.8)[1]]; if (is.na(sam_soft_power)) sam_soft_power <- 8

## 11.5 Construct the co-expression network ------------
# Signed network so only positively co-expressed genes group together; merge close modules.
sub <- ScaleData(sub, features = GetWGCNAGenes(sub), verbose = FALSE)
sub <- ConstructNetwork(sub, soft_power = sam_soft_power, networkType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25,
                        tom_outdir = "hdWGCNA/SAM/TOM", tom_name = "SAM_dev", overwrite_tom = TRUE)
pdf("hdWGCNA/SAM/SAM_dev_module_dendrogram.pdf", width = 10, height = 7)
PlotDendrogram(sub, main = "SAM developmental core hdWGCNA dendrogram"); dev.off()

## 11.6 Module eigengenes and intramodular connectivity ------------
# Harmonize eigengenes across sections (orig.ident) to remove batch structure; kME = hub score.
sub <- ModuleEigengenes(sub, group.by.vars = "orig.ident", assay = "SCT")
sub <- ModuleConnectivity(sub, group.by = "sam_dev_hdwgcna_group", group_name = sam_dev_types)
sub <- ResetModuleNames(sub, new_name = "SAM-M")   # rename modules SAM-M1, SAM-M2, ...

## 11.7 Export modules, hub genes, eigengenes ------------
sam_modules <- GetModules(sub) %>% filter(module != "grey")
sam_hub <- GetHubGenes(sub, n_hubs = 25)
sam_MEs <- GetMEs(sub, harmonized = TRUE)
sam_me_cols <- grep("^SAM-M", colnames(sam_MEs), value = TRUE)
write.csv(sam_hub, "tables/hdWGCNA/SAM_hdWGCNA_hub_genes_2026.7.9.csv", row.names = FALSE)
write.csv(sam_modules, "tables/hdWGCNA/SAM_hdWGCNA_module_assignments_2026.7.9.csv", row.names = FALSE)
sam_module_sizes <- sam_modules %>% dplyr::count(module, name = "n_genes") %>% arrange(desc(n_genes))
write.csv(sam_module_sizes, "tables/hdWGCNA/SAM_hdWGCNA_module_sizes_2026.7.9.csv", row.names = FALSE)

## Plot top hub genes for each module on spatial sections 1, 3, and 4.
sam_hub_gene_col <- if ("gene_name" %in% colnames(sam_hub)) "gene_name" else "gene"
sam_hub_top10 <- sam_hub %>% mutate(gene = .data[[sam_hub_gene_col]]) %>%
  group_by(module) %>% arrange(desc(kME), .by_group = TRUE) %>% dplyr::slice_head(n = 10) %>% ungroup()
write.csv(sam_hub_top10, "tables/hdWGCNA/SAM_hdWGCNA_top10_hub_genes_per_module.csv", row.names = FALSE)
sam_hub_spatial_images <- names(sobj_sam_ann@images)[intersect(c(1, 3, 4), seq_along(names(sobj_sam_ann@images)))]
DefaultAssay(sobj_sam_ann) <- if ("SCT" %in% names(sobj_sam_ann@assays)) "SCT" else DefaultAssay(sobj_sam_ann)
pdf("hdWGCNA/SAM/SAM_hdWGCNA_top10_hub_genes_spatial_images1_3_4.pdf", width = 18, height = 12)
for (module_i in sort(unique(sam_hub_top10$module))) {
  genes_i <- intersect(sam_hub_top10$gene[sam_hub_top10$module == module_i], rownames(sobj_sam_ann))
  if (length(genes_i) == 0) next
  print(SpatialFeaturePlot(sobj_sam_ann, features = genes_i, images = sam_hub_spatial_images, crop = FALSE, ncol = 5,
                           min.cutoff = "q05", max.cutoff = "q95", alpha = c(0.1, 1)) +
          plot_annotation(title = paste0(module_i, " top hub genes on SAM sections 1, 3, and 4")))
}
dev.off()


# Attach metaspot-derived eigengenes using the original spot correspondence.
sub@meta.data[, sam_me_cols] <- sam_MEs[match(rownames(sub@meta.data), rownames(sam_MEs)), sam_me_cols]
## 11.10 GO biological-process enrichment per module ------------
# Normalize GO and gene IDs on both sides; universe = all tested hdWGCNA genes.
go2gene <- read.csv("genome_files/go2gene_poplar.csv", header = FALSE, col.names = c("GO", "gene"), colClasses = "character")
go2term <- read.csv("genome_files/go2term.csv", colClasses = "character", check.names = FALSE)
normalize_go_id <- function(x) ifelse(grepl("^GO:", x), x, paste0("GO:", sprintf("%07d", as.integer(gsub("\\D", "", x)))))
strip_v <- function(g) sub("\\.v5\\.1$", "", as.character(g))
go2gene <- go2gene %>% mutate(GO = normalize_go_id(GO), gene_core = strip_v(gene)) %>% filter(!is.na(GO), !is.na(gene_core), gene_core != "")
go_id_col <- intersect(c("GO", "ID", "id", "go_id"), colnames(go2term))[1]
go_name_col <- intersect(c("name", "Description", "description", "term"), colnames(go2term))[1]
go_class_col <- intersect(c("class", "Ontology", "ontology"), colnames(go2term))[1]
go2term <- go2term %>% mutate(GO = normalize_go_id(.data[[go_id_col]]), name = .data[[go_name_col]], ontology_use = .data[[go_class_col]])
go2term_bp <- go2term %>% filter(ontology_use %in% c("biological_process", "BP"))
go2gene_bp <- go2gene %>% filter(GO %in% go2term_bp$GO) %>% distinct(GO, gene_core)
sam_module_gene_col <- if ("gene" %in% colnames(sam_modules)) "gene" else "gene_name"
sam_go_universe <- strip_v(unique(sam_modules[[sam_module_gene_col]]))
sam_go_overlap <- intersect(sam_go_universe, go2gene_bp$gene_core)
if (length(sam_go_overlap) < 10) stop("Too few hdWGCNA genes map to GO terms. Check gene ID format. Examples module: ", paste(head(sam_go_universe), collapse = ", "), " | GO map: ", paste(head(go2gene_bp$gene_core), collapse = ", "))
sam_go <- lapply(sort(unique(sam_modules$module)), function(m) {
  genes <- strip_v(sam_modules[[sam_module_gene_col]][sam_modules$module == m])
  e <- tryCatch(enricher(genes, TERM2GENE = go2gene_bp[, c("GO", "gene_core")],
                         TERM2NAME = go2term_bp[, c("GO", "name")], universe = sam_go_universe,
                         pvalueCutoff = 0.05, qvalueCutoff = 0.2, maxGSSize = 800),
                error = function(e) NULL)
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  as.data.frame(e) %>% mutate(module = m)
})
sam_go_df <- do.call(rbind, sam_go)
write.csv(sam_go_df, "tables/hdWGCNA/SAM_hdWGCNA_module_GO_2026.7.9.csv", row.names = FALSE)

## Plot hdWGCNA module GO enrichment, one module per page.
if (!is.null(sam_go_df) && nrow(sam_go_df) > 0) {
  sam_go_plot_df <- sam_go_df %>% group_by(module) %>% arrange(p.adjust, .by_group = TRUE) %>% slice_head(n = 12) %>% ungroup() %>%
    mutate(Description = stringr::str_trunc(Description, 70),
           GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x))),
           neg_log10_padj = -log10(pmax(p.adjust, 1e-300)))
  sam_go_color_limit <- quantile(sam_go_plot_df$neg_log10_padj, 0.95, na.rm = TRUE)
  p_sam_go_all <- ggplot(sam_go_plot_df, aes(GeneRatio_num, forcats::fct_reorder(Description, GeneRatio_num))) +
    geom_point(aes(size = Count, color = neg_log10_padj)) + facet_wrap(~module, scales = "free_y", ncol = 3) +
    scale_color_gradientn(colors = c("grey85", "#FCAE91", "#FB6A4A", "#CB181D"), limits = c(0, sam_go_color_limit), oob = scales::squish) +
    theme_classic() + labs(title = "SAM hdWGCNA module GO biological processes", x = "Gene ratio", y = NULL, color = "-log10(adj. P)", size = "Gene count")
  ggsave("hdWGCNA/SAM/SAM_hdWGCNA_module_GO_dotplot_all_modules.pdf", p_sam_go_all, width = 14, height = 12)
  pdf("hdWGCNA/SAM/SAM_hdWGCNA_module_GO_dotplot_per_module.pdf", width = 8, height = 6)
  for (module_i in sort(unique(sam_go_plot_df$module))) {
    go_i <- sam_go_plot_df %>% filter(module == module_i) %>% arrange(p.adjust, desc(GeneRatio_num)) %>% mutate(Description = forcats::fct_reorder(Description, GeneRatio_num))
    p_go_i <- ggplot(go_i, aes(GeneRatio_num, Description)) + geom_point(aes(size = Count, color = neg_log10_padj)) +
      scale_color_gradientn(colors = c("grey85", "#FCAE91", "#FB6A4A", "#CB181D"), limits = c(0, sam_go_color_limit), oob = scales::squish) +
      theme_classic() + labs(title = paste0(module_i, " GO biological process"), x = "Gene ratio", y = NULL, color = "-log10(adj. P)", size = "Gene count")
    print(p_go_i)
  }
  dev.off()
}


for (nm in sam_me_cols) {
  target <- paste0("ME_", nm)
  sobj_sam_ann[[target]] <- NA_real_
  sobj_sam_ann@meta.data[colnames(sub), target] <- sub@meta.data[[nm]]
}
qs::qsave(sub, "saved_obj/sobj_sam_dev_hdwgcna.qs")
qs::qsave(sobj_sam_ann, "saved_obj/sobj_sam_res0.4_final_v1_2026.7.24.qs")

  },
  trichome = function() {
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

# Figure contract
# Conclusion: independently supported poplar trichome markers define initiation,
# pan-trichome, and developing programs that localize to leaf primordia and change
# along the existing meristem-to-primordium trajectory.
# Analysis scope: reuse the final spatial object, scRNA-seq reference, RCTD weights,
# trajectory, NMF factors, and hdWGCNA eigengenes. Recompute marker-derived scores
# and all correlations that depend on those scores.

set.seed(12345)

workspace_dir <- file.path(ATLAS_WORK_ROOT)
supplement_file <- file.path(ATLAS_WORK_ROOT, "reference", "Giabardo_2026_Supplemental_Tables.xlsx")
spatial_file <- file.path(workspace_dir, "saved_obj", "sobj_sam_res0.4_final_v1_2026.7.24.qs")
scrna_file <- file.path(workspace_dir, "saved_obj", "scRNA", "sam_cleaned_downsampled_12k_SCT_for_spatial_transfer.qs")
rctd_file <- file.path(workspace_dir, "deconvolution", "RCTD", "tables", "SAM_RCTD_multi_celltype_weights.csv")

revision_tag <- "2026-09-14"
out_dir <- file.path(workspace_dir, "figures", "Shoot_apex", paste0("Figure4_trichome_revision_", revision_tag))
table_dir <- file.path(workspace_dir, "tables", paste0("Figure4_trichome_revision_", revision_tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(supplement_file), file.exists(spatial_file), file.exists(scrna_file), file.exists(rctd_file))

palette_contract <- c(
  initial = "#2C7FB8",
  pan = "#7B3294",
  developing = "#D95F0E",
  core = "#1B7837",
  neutral = "#737373"
)

theme_pub <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.6),
      strip.text = element_text(size = base_size - 0.1, face = "bold"),
      plot.title = element_text(size = base_size + 0.8, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank()
    )
}

theme_set(theme_pub())

save_pub <- function(plot, stem, width_mm = 183, height_mm = 120, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 300, background = "white")
  print(plot)
  dev.off()
}

resolve_genes <- function(ids, universe) {
  vapply(ids, function(id) {
    hit <- intersect(c(id, paste0(id, ".v5.1")), universe)
    if (length(hit)) hit[[1]] else NA_character_
  }, character(1))
}

zscore_gene_mean <- function(object, genes, assay = "SCT") {
  genes <- na.omit(resolve_genes(genes, rownames(object[[assay]])))
  if (!length(genes)) return(rep(NA_real_, ncol(object)))
  x <- as.matrix(GetAssayData(object, assay = assay, layer = "data")[genes, , drop = FALSE])
  z <- t(scale(t(x)))
  z[!is.finite(z)] <- 0
  colMeans(z)
}

# The aliases below follow the companion paper and its supplementary tables.
# They deliberately replace the earlier BLAST-only MYB38 assignments.
marker_tbl <- tribble(
  ~gene_id,                ~gene_name, ~program,     ~evidence, ~include_in_core,
  "PtXaAlbH.10G130200",    "MYB38",   "Initial",   "Companion Supplementary Table 2; glabrous-mutant DE", TRUE,
  "PtXaTreH.10G136900",    "MYB38",   "Initial",   "Reciprocal 717 allele; glabrous-mutant DE", TRUE,
  "PtXaTreH.10G194200",    "SMR1",    "Pan",       "Companion Supplementary Table 3; all trichome subclusters", TRUE,
  "PtXaAlbH.10G187700",    "SMR1",    "Pan",       "Second allele named SMR1 in companion paper; glabrous-mutant DE", TRUE,
  "PtXaAlbH.10G166200",    "MATL16",  "Developing","Companion paper; developing trichomes; promoter-supported gene pair", TRUE,
  "PtXaTreH.10G173100",    "MATL16",  "Developing","Companion Supplementary Table 3; promoter validated", TRUE,
  "PtXaAlbH.19G048000",    "MYB59",   "Developing","Companion paper; developing trichome subcluster", FALSE,
  "PtXaTreH.19G107500",    "LAC3",    "Developing","Companion paper; also detected in xylem fiber", FALSE
)

bulk_de <- read_excel(supplement_file, sheet = "Supp. Table 7 - Bulk DE_update", skip = 2) |>
  transmute(
    gene_id = .data[["Poplar 717 Gene ID"]],
    bulk_annotation = .data[["Functional Annotation"]],
    best_arabidopsis_hit = .data[["Best Aarabidopsis BLAST hit"]],
    ko_vs_control_log2fc = as.numeric(.data[["KO vs Control log2FC"]]),
    ko_vs_control_padj = as.numeric(.data[["KO vs Control Adjusted p-value"]])
  )

marker_tbl <- marker_tbl |>
  left_join(bulk_de, by = "gene_id")

scrna <- qread(scrna_file)
spatial <- qread(spatial_file)
DefaultAssay(scrna) <- "RNA"
DefaultAssay(spatial) <- "SCT"

marker_tbl <- marker_tbl |>
  mutate(
    scrna_feature = resolve_genes(gene_id, rownames(scrna[["RNA"]])),
    spatial_feature = resolve_genes(gene_id, rownames(spatial[["SCT"]])),
    present_scrna = !is.na(scrna_feature),
    present_spatial = !is.na(spatial_feature)
  )

write.csv(marker_tbl, file.path(table_dir, "companion_marker_gene_audit.csv"), row.names = FALSE)

if (sum(marker_tbl$present_spatial & marker_tbl$include_in_core) < 4) {
  stop("Fewer than four core companion markers were found in the spatial object.")
}

score_sets <- list(
  "Core identity" = marker_tbl |> filter(include_in_core) |> pull(gene_id),
  "Initial (MYB38)" = marker_tbl |> filter(program == "Initial") |> pull(gene_id),
  "Pan-trichome (SMR1)" = marker_tbl |> filter(program == "Pan") |> pull(gene_id),
  "Developing (MATL16/MYB59)" = marker_tbl |> filter(program == "Developing", gene_name != "LAC3") |> pull(gene_id)
)

score_columns <- c(
  "Core identity" = "trichome_core_score",
  "Initial (MYB38)" = "trichome_initial_score",
  "Pan-trichome (SMR1)" = "trichome_pan_score",
  "Developing (MATL16/MYB59)" = "trichome_developing_score"
)

for (nm in names(score_sets)) {
  spatial[[score_columns[[nm]]]] <- zscore_gene_mean(spatial, score_sets[[nm]], assay = "SCT")
}

# Import independent RCTD evidence for the two single-cell clusters identified by
# the companion study as trichome initials (cluster 14) and developing trichomes
# (cluster 39). These are evidence layers, not ingredients in the marker score.
rctd <- read.csv(rctd_file, check.names = FALSE)
rownames(rctd) <- rctd$spot
rctd$spot <- NULL
common_spots <- intersect(colnames(spatial), rownames(rctd))
spatial$rctd_trichome_initial <- NA_real_
spatial$rctd_trichome_developing <- NA_real_
spatial$rctd_trichome_initial[match(common_spots, colnames(spatial))] <- rctd[common_spots, "14: Epidermis"]
spatial$rctd_trichome_developing[match(common_spots, colnames(spatial))] <- rctd[common_spots, "39: Epidermis"]

object_target <- file.path(workspace_dir, "saved_obj", paste0("sobj_shoot_apex_trichome_scores_companion_", revision_tag, ".qs"))
object_temp <- tempfile(pattern = "trichome_scores_", tmpdir = tempdir(), fileext = ".qs")
qsave(spatial, object_temp)
stopifnot(file.copy(object_temp, object_target, overwrite = TRUE))

# scRNA-seq expression audit across all reference clusters.
scrna_genes <- na.omit(marker_tbl$scrna_feature)
scrna_mat <- GetAssayData(scrna, assay = "RNA", layer = "data")[scrna_genes, , drop = FALSE]
scrna_cluster <- as.character(Idents(scrna))
scrna_levels <- names(sort(table(scrna_cluster), decreasing = TRUE))
scrna_summary <- bind_rows(lapply(scrna_levels, function(cluster_i) {
  keep <- scrna_cluster == cluster_i
  tibble(
    scrna_feature = scrna_genes,
    cluster = cluster_i,
    average_expression = Matrix::rowMeans(scrna_mat[, keep, drop = FALSE]),
    percent_expressing = 100 * Matrix::rowMeans(scrna_mat[, keep, drop = FALSE] > 0),
    n_cells = sum(keep)
  )
})) |>
  left_join(marker_tbl |> select(scrna_feature, gene_id, gene_name, program), by = "scrna_feature")

write.csv(scrna_summary, file.path(table_dir, "companion_markers_scrna_expression_by_cluster.csv"), row.names = FALSE)

focus_clusters <- c(
  "14: Epidermis", "39: Epidermis", "19: Epidermis", "18: Epidermis",
  "4: Epidermis", "22: Dividing epidermis", "24: Vessel elements"
)

gene_axis <- marker_tbl |>
  filter(present_scrna) |>
  mutate(label = paste0(gene_name, "\n", gene_id), label = factor(label, levels = rev(label))) |>
  select(scrna_feature, label)

dot_df <- scrna_summary |>
  filter(cluster %in% focus_clusters) |>
  inner_join(gene_axis, by = "scrna_feature") |>
  group_by(scrna_feature) |>
  mutate(relative_expression = if (max(average_expression) > 0) average_expression / max(average_expression) else 0) |>
  ungroup() |>
  mutate(cluster = factor(cluster, levels = focus_clusters))

write.csv(dot_df, file.path(table_dir, "companion_markers_scrna_dotplot_source_data.csv"), row.names = FALSE)

p_fig4c <- ggplot(dot_df, aes(cluster, label)) +
  geom_point(aes(size = percent_expressing, colour = relative_expression), stroke = 0.2) +
  scale_size_continuous(range = c(0.2, 5.5), breaks = c(10, 25, 50, 75), limits = c(0, 100)) +
  scale_colour_gradient(low = "#F2F2F2", high = "#2166AC", limits = c(0, 1)) +
  labs(
    title = "Poplar single-cell atlas identifies trichome-associated states",
    subtitle = "Clusters 14 and 39 represent trichome initials and developing trichomes, respectively",
    x = NULL, y = NULL, size = "Expressing cells (%)", colour = "Relative mean"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 38, hjust = 1), axis.text.y = element_text(face = "italic"), legend.position = "bottom")

save_pub(p_fig4c, file.path(out_dir, "Fig4C_companion_scRNA_trichome_markers"), width_mm = 183, height_mm = 112)

# Spatial score maps and independent deconvolution maps.
representative_section <- "sam_A_s1"
if (!representative_section %in% Images(spatial)) representative_section <- Images(spatial)[1]
spatial_features <- c("trichome_core_score", "rctd_trichome_initial", "rctd_trichome_developing")
spatial_titles <- c(
  "Trichome identity score",
  "scRNA cluster 14\ntrichome initials",
  "scRNA cluster 39\ndeveloping trichomes"
)
spatial_legend_titles <- c("Score (gene z-score)", "RCTD weight", "RCTD weight")
spatial_plots <- SpatialFeaturePlot(
  spatial,
  features = spatial_features,
  images = representative_section,
  crop = FALSE,
  combine = FALSE,
  pt.size.factor = 1.75,
  min.cutoff = "q02",
  max.cutoff = "q98"
)
spatial_plots <- Map(function(p, ttl, legend_title) {
  p +
    labs(title = ttl, fill = legend_title, colour = legend_title) +
    guides(
      fill = guide_colorbar(title.position = "top", barwidth = grid::unit(24, "mm")),
      colour = guide_colorbar(title.position = "top", barwidth = grid::unit(24, "mm"))
    ) +
    theme_pub() +
    theme(
      axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      legend.position = "bottom", legend.title = element_text(size = 6),
      plot.margin = margin(2, 5, 2, 5)
    )
}, spatial_plots, spatial_titles, spatial_legend_titles)
p_fig4d <- wrap_plots(spatial_plots, nrow = 1) +
  plot_annotation(
    title = "Validated trichome markers and single-cell states converge spatially",
    subtitle = "Representative shoot-apex section",
    tag_levels = "a"
  ) & theme(plot.tag = element_text(size = 8, face = "bold"))

save_pub(p_fig4d, file.path(out_dir, "Fig4D_spatial_marker_score_and_deconvolution"), width_mm = 183, height_mm = 78)

# Marker scores along the existing trajectory. The trajectory itself is unchanged.
meta <- spatial@meta.data |>
  rownames_to_column("spot")
subcluster_levels <- c(paste0("Leaf primordium", 1:6), "Meristem")
subcluster_colors <- setNames(scales::hue_pal()(7), subcluster_levels)
write.csv(data.frame(subcluster = subcluster_levels, color = unname(subcluster_colors)),
          file.path(table_dir, "Figure4a_subcluster_palette.csv"), row.names = FALSE)
upward <- meta |>
  filter(is.finite(monocle3_pseudotime_upward_leaf)) |>
  mutate(
    subcluster = factor(gsub("\\s+", " ", trimws(monocle3_subcluster_upward_leaf)), levels = subcluster_levels)
  )

score_long <- upward |>
  select(spot, pseudotime = monocle3_pseudotime_upward_leaf, subcluster, all_of(unname(score_columns))) |>
  pivot_longer(all_of(unname(score_columns)), names_to = "score_column", values_to = "score") |>
  mutate(
    signature = factor(names(score_columns)[match(score_column, score_columns)], levels = names(score_sets))
  )

score_deciles <- score_long |>
  group_by(signature) |>
  mutate(decile = ntile(pseudotime, 10)) |>
  group_by(signature, decile) |>
  summarize(
    pseudotime_mid = median(pseudotime),
    mean_score = mean(score),
    se_score = sd(score) / sqrt(n()),
    n_spots = n(),
    .groups = "drop"
  )
write.csv(score_deciles, file.path(table_dir, "trichome_stage_scores_by_pseudotime_decile.csv"), row.names = FALSE)

stopifnot(!anyNA(upward$subcluster))
write.csv(score_long |> filter(signature == "Core identity"),
          file.path(table_dir, "core_score_pseudotime_by_subcluster.csv"), row.names = FALSE)

p_fig4e <- ggplot(score_long |> filter(signature == "Core identity"), aes(pseudotime, score)) +
  geom_point(aes(colour = subcluster), size = 0.55, alpha = 0.75) +
  geom_smooth(method = "loess", se = TRUE, span = 0.7, colour = "black", fill = "grey75", linewidth = 0.8) +
  scale_colour_manual(values = subcluster_colors, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
  labs(
    title = "Poplar trichome identity changes along the existing developmental trajectory",
    subtitle = "Core score uses poplar MYB38, SMR1, and MATL16 alleles (Giabardo et al., 2026)",
    x = "Meristem-to-leaf-primordium pseudotime", y = "Trichome core score (mean gene z-score)", colour = "Subcluster"
  ) +
  theme_pub() + theme(legend.position = "right")

save_pub(p_fig4e, file.path(out_dir, "Fig4E_core_trichome_score_along_pseudotime"), width_mm = 128, height_mm = 90)

p_stage <- ggplot(score_deciles, aes(pseudotime_mid, mean_score, colour = signature, fill = signature)) +
  geom_ribbon(aes(ymin = mean_score - se_score, ymax = mean_score + se_score), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(
    "Core identity" = palette_contract[["core"]],
    "Initial (MYB38)" = palette_contract[["initial"]],
    "Pan-trichome (SMR1)" = palette_contract[["pan"]],
    "Developing (MATL16/MYB59)" = palette_contract[["developing"]]
  )) +
  scale_fill_manual(values = c(
    "Core identity" = palette_contract[["core"]],
    "Initial (MYB38)" = palette_contract[["initial"]],
    "Pan-trichome (SMR1)" = palette_contract[["pan"]],
    "Developing (MATL16/MYB59)" = palette_contract[["developing"]]
  )) +
  labs(
    title = "Stage-resolved trichome signatures along pseudotime",
    x = "Meristem-to-leaf-primordium pseudotime", y = "Mean score", colour = NULL, fill = NULL
  ) +
  theme_pub() + theme(legend.position = "bottom")

save_pub(p_stage, file.path(out_dir, "FigS_stage_resolved_trichome_scores_along_pseudotime"), width_mm = 128, height_mm = 88)

# Recompute module-trait correlations. The hdWGCNA modules are retained because
# their construction does not depend on the old trichome marker score.
module_cols <- grep("^ME_SAM-M", colnames(meta), value = TRUE)
trait_cols <- c(
  trichome_core_score = "Core trichome score",
  trichome_initial_score = "MYB38 initiation score",
  trichome_pan_score = "SMR1 pan-trichome score",
  trichome_developing_score = "MATL16/MYB59 developing score",
  rctd_trichome_initial = "Predicted scRNA-seq cluster 14 contribution (trichome initials)",
  rctd_trichome_developing = "Predicted scRNA-seq cluster 39 contribution (developing trichomes)",
  monocle3_pseudotime_upward_leaf = "Meristem-to-primordium pseudotime"
)

dev_meta <- meta |>
  filter(celltypes %in% c("Meristem", "Leaf primordium"))

cor_df <- bind_rows(lapply(module_cols, function(module_i) {
  bind_rows(lapply(names(trait_cols), function(trait_i) {
    ok <- is.finite(dev_meta[[module_i]]) & is.finite(dev_meta[[trait_i]])
    if (sum(ok) < 10 || sd(dev_meta[[module_i]][ok]) == 0 || sd(dev_meta[[trait_i]][ok]) == 0) {
      return(tibble(module = sub("^ME_", "", module_i), trait = trait_cols[[trait_i]], correlation = NA_real_, p_value = NA_real_, n = sum(ok)))
    }
    ct <- suppressWarnings(cor.test(dev_meta[[module_i]][ok], dev_meta[[trait_i]][ok], method = "pearson"))
    tibble(module = sub("^ME_", "", module_i), trait = trait_cols[[trait_i]], correlation = unname(ct$estimate), p_value = ct$p.value, n = sum(ok))
  }))
})) |>
  group_by(trait) |>
  mutate(p_adjust = p.adjust(p_value, method = "BH")) |>
  ungroup() |>
  mutate(significance = case_when(
    p_adjust < 0.001 ~ "***",
    p_adjust < 0.01 ~ "**",
    p_adjust < 0.05 ~ "*",
    TRUE ~ ""
  ))

write.csv(cor_df, file.path(table_dir, "updated_hdwgcna_module_trait_correlations.csv"), row.names = FALSE)

module_order <- cor_df |>
  filter(trait == "Core trichome score") |>
  arrange(correlation) |>
  pull(module)
trait_order <- unname(trait_cols)
p_module <- cor_df |>
  mutate(module = factor(module, levels = module_order), trait = factor(trait, levels = rev(trait_order))) |>
  ggplot(aes(module, trait, fill = correlation)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = significance), size = 2.2, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
  labs(
    title = "Co-expression modules associated with poplar trichome markers",
    subtitle = "Pearson correlations across meristem and leaf-primordium spots; stars show BH-adjusted P values",
    x = NULL, y = NULL, fill = "Pearson r"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

save_pub(p_module, file.path(out_dir, "FigS_full_module_trait_correlations"), width_mm = 220, height_mm = 105)
main_traits <- unname(trait_cols[c("trichome_core_score", "rctd_trichome_initial",
                                  "rctd_trichome_developing", "monocle3_pseudotime_upward_leaf")])
main_cor_df <- cor_df |> filter(trait %in% main_traits) |>
  mutate(module = factor(module, levels = module_order),
         trait = factor(trait, levels = rev(main_traits)))
write.csv(main_cor_df, file.path(table_dir, "main_figure_module_trait_correlations.csv"), row.names = FALSE)
p_main_module <- p_module %+% main_cor_df +
  scale_y_discrete(labels = function(x) gsub(" contribution ", " contribution\n", x, fixed = TRUE)) +
  labs(title = "Co-expression modules associated with trichome identity and development",
       subtitle = "Predicted scRNA-seq contributions estimated by RCTD; stars show BH-adjusted P values")
save_pub(p_main_module, file.path(out_dir, "Fig4F_updated_module_trait_correlations"), width_mm = 210, height_mm = 76)

# Quantitative audit: domain enrichment, agreement with RCTD, and score-trajectory correlations.
domain_summary <- meta |>
  group_by(celltypes) |>
  summarize(
    n_spots = n(),
    across(all_of(c(unname(score_columns), "rctd_trichome_initial", "rctd_trichome_developing")),
           list(mean = ~mean(.x, na.rm = TRUE), median = ~median(.x, na.rm = TRUE))),
    .groups = "drop"
  )
write.csv(domain_summary, file.path(table_dir, "trichome_scores_and_rctd_by_spatial_domain.csv"), row.names = FALSE)

metric_pairs <- tribble(
  ~x, ~y, ~comparison,
  "trichome_core_score", "rctd_trichome_initial", "Core score vs RCTD cluster 14",
  "trichome_core_score", "rctd_trichome_developing", "Core score vs RCTD cluster 39",
  "trichome_initial_score", "rctd_trichome_initial", "MYB38 score vs RCTD cluster 14",
  "trichome_developing_score", "rctd_trichome_developing", "Developing score vs RCTD cluster 39",
  "trichome_core_score", "monocle3_pseudotime_upward_leaf", "Core score vs upward pseudotime"
)
if ("trichome_score" %in% colnames(meta)) {
  metric_pairs <- bind_rows(
    metric_pairs,
    tibble(
      x = "trichome_core_score", y = "trichome_score",
      comparison = "Companion core score vs previous manuscript score"
    )
  )
}

agreement <- metric_pairs |>
  rowwise() |>
  mutate(
    n = sum(is.finite(meta[[x]]) & is.finite(meta[[y]])),
    spearman_rho = suppressWarnings(cor(meta[[x]], meta[[y]], method = "spearman", use = "pairwise.complete.obs")),
    pearson_r = suppressWarnings(cor(meta[[x]], meta[[y]], method = "pearson", use = "pairwise.complete.obs"))
  ) |>
  ungroup()
write.csv(agreement, file.path(table_dir, "trichome_score_agreement_statistics.csv"), row.names = FALSE)

session <- capture.output(sessionInfo())
writeLines(session, file.path(table_dir, "sessionInfo.txt"))

message("Completed targeted Figure 4 trichome reanalysis.")
message("Figure outputs: ", out_dir)
message("Source-data tables: ", table_dir)

saveRDS(list(pA=p_fig4c, pB=spatial_plots, pE=p_fig4e, pF=p_main_module,
 pStage=p_stage, pFull=p_module, upward=upward, marker_tbl=marker_tbl,
 subcluster_colors=subcluster_colors, meta=meta), "saved_obj/figure4_plot_inputs.rds")

  }
)
atlas_dispatch(stages)
