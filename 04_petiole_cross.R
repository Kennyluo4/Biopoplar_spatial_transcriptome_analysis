# Submission edition: run from this directory or set POPLAR_CODE_ROOT.
.code_root <- Sys.getenv("POPLAR_CODE_ROOT", unset = "")
if (!nzchar(.code_root)) {
  .script <- grep("^--file=", commandArgs(), value = TRUE)
  .code_root <- if (length(.script)) dirname(normalizePath(sub("^--file=", "", .script[1]))) else getwd()
}
source(file.path(.code_root, "R", "submission_setup.R"))

stages <- list(
  preprocess = function() {


# Improved spatial Seurat pipeline for poplar petiole cross-sections
# - Automatically split disconnected sections with DBSCAN
# - Merge all resulting sections
# - Run SCT + Harmony integration + clustering
# - Generate marker-based provisional cell type annotation
# - Save objects, plots, and marker tables
#
##linux R version: R/4.3.2-gfbf-2023a

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(harmony)
  library(scCustomize)
  library(qs)
  library(ComplexHeatmap)
  library(dbscan)
  library(psych)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(circlize)
})
dir.create("QC", showWarnings = FALSE, recursive = TRUE)
dir.create("clustering", showWarnings = FALSE, recursive = TRUE)
dir.create("marker", showWarnings = FALSE, recursive = TRUE)
dir.create("saved_obj", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

# ----------------------------- #
# 0. Helper functions ###########
# ----------------------------- #

load_slice <- function(data_dir, slice_name) {
  obj <- Load10X_Spatial(data.dir = data_dir, slice = slice_name)
  obj$orig.ident <- slice_name
  obj$slice_id <- slice_name
  obj
}

detect_sections_dbscan <- function(seu, eps = 500, minPts = 10, xcol = "x", ycol = "y") {
  coords <- GetTissueCoordinates(seu)
  stopifnot(all(c(xcol, ycol) %in% colnames(coords)))
  mat <- as.matrix(coords[, c(xcol, ycol)])
  db <- dbscan::dbscan(mat, eps = eps, minPts = minPts)
  section_auto <- db$cluster
  names(section_auto) <- rownames(coords)
  section_auto
}

# split a Seurat spatial object into one object per disconnected tissue section
# objects are ordered top-to-bottom then left-to-right based on centroids
split_spatial_sections <- function(
  seu,
  base_name,
  eps = 500,
  minPts = 10,
  xcol = "x",
  ycol = "y",
  keep_noise = FALSE
) {
  coords <- GetTissueCoordinates(seu)
  sec <- detect_sections_dbscan(seu, eps = eps, minPts = minPts, xcol = xcol, ycol = ycol)

  valid_clusters <- sort(unique(sec))
  if (!keep_noise) {
    valid_clusters <- valid_clusters[valid_clusters != 0]
  }

  if (length(valid_clusters) == 0) {
    stop("No DBSCAN sections detected for ", base_name, ". Try increasing eps or decreasing minPts.")
  }

  centroids <- lapply(valid_clusters, function(k) {
    idx <- names(sec)[sec == k]
    tibble(
      cluster = k,
      x_center = median(coords[idx, xcol]),
      y_center = median(coords[idx, ycol]),
      n_spots = length(idx)
    )
  }) %>% bind_rows()

  # order top-to-bottom, then left-to-right
  centroids <- centroids %>%
    arrange(y_center, x_center) %>%
    mutate(section_rank = row_number(),
           new_name = paste0(base_name, "_s", section_rank))

  out <- vector("list", nrow(centroids))
  names(out) <- centroids$new_name

  original_image_name <- names(seu@images)[1]

  for (i in seq_len(nrow(centroids))) {
    cl <- centroids$cluster[i]
    new_name <- centroids$new_name[i]
    keep_cells <- names(sec)[sec == cl]
    sub_obj <- subset(seu, cells = keep_cells)

    # rename image slot to the new section name
    names(sub_obj@images) <- new_name
    sub_obj$orig.ident <- new_name
    sub_obj$slice_id <- new_name
    sub_obj$parent_library <- base_name
    sub_obj$section_rank <- centroids$section_rank[i]
    sub_obj$section_auto <- paste0("section_", centroids$section_rank[i])
    out[[new_name]] <- sub_obj
  }

  attr(out, "centroids") <- centroids
  attr(out, "original_image_name") <- original_image_name
  out
}

plot_dbscan_check <- function(seu, title = NULL, pt.size.factor = 1.8) {
  coords <- GetTissueCoordinates(seu)
  seu$section_auto <- detect_sections_dbscan(seu)
  p1 <- SpatialDimPlot(seu, group.by = "section_auto", pt.size.factor = pt.size.factor) +
    ggtitle(ifelse(is.null(title), seu$orig.ident[1], title))
  p2 <- ggplot(coords, aes(x = x, y = y, color = factor(seu$section_auto[rownames(coords)]))) +
    geom_point(size = 1) +
    scale_y_reverse() +
    coord_equal() +
    theme_classic() +
    labs(color = "DBSCAN", title = "Coordinate view")
  p1 | p2
}

process_spatial_integration <- function(sobj, resolution = 0.4, dims = 1:30) {
  sobj <- PercentageFeatureSet(sobj, pattern = "^ChrMt", col.name = "percent_mt")
  sobj <- PercentageFeatureSet(sobj, pattern = "^ChrPt", col.name = "percent_chlp")

  sobj <- SCTransform(sobj, assay = "Spatial", verbose = FALSE)
  sobj <- RunPCA(sobj, assay = "SCT", verbose = FALSE)
  sobj <- RunHarmony(
    object = sobj,
    group.by.vars = "orig.ident",
    reduction.save = "harmony",
    plot_convergence = FALSE
  )
  sobj <- RunUMAP(sobj, reduction = "harmony", dims = dims)
  sobj <- FindNeighbors(sobj, reduction = "harmony", dims = dims)
  sobj <- FindClusters(sobj, resolution = resolution)
  sobj
}

plot_cluster_proportions <- function(sobj, order_vec, cluster_col = NULL, out_pdf) {
  if (is.null(cluster_col)) {
    plot_data <- sobj@meta.data %>%
      mutate(cluster = as.character(Idents(sobj)))
  } else {
    plot_data <- sobj@meta.data %>%
      mutate(cluster = as.character(.data[[cluster_col]]))
  }
  
  plot_data <- plot_data %>%
    group_by(orig.ident, cluster) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(orig.ident) %>%
    mutate(proportion = n_cells / sum(n_cells))
  
  plot_data$orig.ident <- factor(plot_data$orig.ident, levels = order_vec)
  
  p <- ggplot(plot_data, aes(x = orig.ident, y = proportion, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = cluster), position = position_stack(vjust = 0.5),
              size = 3, color = "white") +
    theme_classic() +
    labs(x = NULL, y = "Proportion", fill = ifelse(is.null(cluster_col), "Ident", cluster_col)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(out_pdf, p, width = 11, height = 6)
}

make_marker_heatmap <- function(sobj, marker_df, out_pdf, group_col = "celltypes") {
  Idents(sobj) <- group_col
  
  df_avg <- AverageExpression(sobj, assays = "SCT", layer = "data")$SCT %>% as.matrix()
  keep_var <- apply(df_avg, 1, var) > 0
  z_mtx <- t(scale(t(df_avg[keep_var, , drop = FALSE])))
  z_mtx[!is.finite(z_mtx)] <- 0
  
  plot_data <- marker_df %>%
    filter(GeneID.v5 %in% rownames(z_mtx)) %>%
    mutate(
      row_label = ifelse(
        is.na(Name) | Name == "",
        GeneID.v5,
        paste0(Name, " (", GeneID.v5, ")")
      )
    )
  
  ht_matrix <- z_mtx[plot_data$GeneID.v5, , drop = FALSE]
  rownames(ht_matrix) <- plot_data$row_label
  
  pdf(out_pdf, height = 18, width = 10)
  ht <- Heatmap(
    ht_matrix,
    name = "Z-score",
    row_split = plot_data$celltype,
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 8),
    row_gap = grid::unit(3, "mm"),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 5.5),
    column_names_gp = grid::gpar(fontsize = 9),
    heatmap_legend_param = list(direction = "horizontal")
  )
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()
}

score_and_assign_celltypes <- function(sobj, marker_df, assay = "SCT", cluster_col = "seurat_clusters") {
  DefaultAssay(sobj) <- assay
  
  # remove old score columns if present
  md <- sobj@meta.data
  md <- md[, !grepl("^score_|^tmpScore_", colnames(md)), drop = FALSE]
  sobj@meta.data <- md
  
  marker_sets <- marker_df %>%
    distinct(celltype, GeneID.v5) %>%
    filter(GeneID.v5 %in% rownames(sobj)) %>%
    group_by(celltype) %>%
    summarise(genes = list(unique(GeneID.v5)), .groups = "drop")
  
  # make safe unique names
  score_names <- make.unique(paste0("score_", make.names(marker_sets$celltype)))
  
  sobj <- AddModuleScore(
    object = sobj,
    features = marker_sets$genes,
    name = "tmpScore_",
    assay = assay,
    search = FALSE
  )
  
  added_cols <- paste0("tmpScore_", seq_len(nrow(marker_sets)))
  names(sobj@meta.data)[match(added_cols, names(sobj@meta.data))] <- score_names
  
  cluster_scores <- sobj@meta.data %>%
    mutate(cluster_label = as.character(.data[[cluster_col]])) %>%
    group_by(cluster_label) %>%
    summarise(across(all_of(score_names), mean), .groups = "drop")
  
  score_long <- cluster_scores %>%
    tidyr::pivot_longer(cols = all_of(score_names), names_to = "score_name", values_to = "mean_score") %>%
    mutate(celltype = marker_sets$celltype[match(score_name, score_names)]) %>%
    group_by(cluster_label) %>%
    arrange(desc(mean_score), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  top_assign <- score_long %>%
    filter(rank == 1) %>%
    select(cluster_label, predicted_celltype = celltype, predicted_score = mean_score)
  
  sobj$predicted_celltype <- top_assign$predicted_celltype[
    match(as.character(sobj@meta.data[[cluster_col]]), top_assign$cluster_label)
  ]
  
  list(
    object = sobj,
    cluster_scores = cluster_scores,
    top_assign = top_assign,
    score_long = score_long
  )
}

plot_annotation_qc <- function(sobj, out_pdf) {
  p1 <- DimPlot(sobj, reduction = "umap", group.by = "orig.ident", alpha = 0.7) + ggtitle("By section")
  p2 <- DimPlot(sobj, reduction = "umap", label = TRUE) + ggtitle("By cluster")
  p3 <- SpatialDimPlot(sobj, crop = FALSE, ncol = 4)
  p4 <- SpatialDimPlot(sobj, group.by = "predicted_celltype", crop = FALSE, ncol = 4)
  pdf(out_pdf, width = 16, height = 10)
  print(p1 + p2)
  print(p3)
  print(p4)
  dev.off()
}

plot_known_markers <- function(sobj, marker_genes, out_pdf, ncol = 4, pt.size.factor = 1.8) {
  marker_genes <- intersect(marker_genes, rownames(sobj))
  if (length(marker_genes) == 0) {
    message("No valid marker genes found for: ", out_pdf)
    return(invisible(NULL))
  }
  pdf(out_pdf, width = 14, height = max(8, ceiling(length(marker_genes) / ncol) * 3))
  print(
    SpatialFeaturePlot(
      object = sobj,
      features = marker_genes,
      alpha = c(0.1, 1),
      ncol = ncol,
      pt.size.factor = pt.size.factor,
      crop = FALSE
    )
  )
  dev.off()
}

rename_clusters_from_table <- function(sobj, mapping_df) {
  mapping_df <- mapping_df %>%
    mutate(seurat_clusters = as.character(cluster_label))
  cl_to_type <- setNames(mapping_df$celltype_final, mapping_df$cluster_label)
  sobj <- RenameIdents(sobj, !!!cl_to_type)
  sobj$celltypes <- Idents(sobj)
  sobj
}


# --------------------------------------- #
# Reviewed common helper overrides ####
# --------------------------------------- #
# This second source call intentionally overrides older local helper definitions
# with shared, reviewed versions. Tissue-specific analysis blocks below are kept.
source(file.path(ATLAS_CODE_ROOT, "R", "common_helpers.R"), local = TRUE)
make_output_dirs()

# Plot layout knobs: adjust here without editing every SpatialDimPlot call.
plot_cfg <- list(
  resolution_width = 14, resolution_height = 12,
  final_width = 16, final_height = 10,
  marker_width = 16, marker_height = 12,
  ncol = 4, pt.size.factor = 1.8, crop = FALSE
)

# ----------------------------- #
# 1. Load raw libraries ######### 
# ----------------------------- #

l1c   <- load_slice("data/petiole1_A_sub1/",   "petiole_l1_cross")
l2c   <- load_slice("data/petiole1_B_sub2/",   "petiole_l2_cross")
l4c_1 <- load_slice("data/petiole1_C_sub2/",   "petiole_l4_cross1")
l4c_2 <- load_slice("data/petiole3_D_L4cross/","petiole_l4_cross2")
l5c_1 <- load_slice("data/petiole1_D_sub2/",   "petiole_l5_cross1")
l5c_2 <- load_slice("data/petiole3_D_L5cross/","petiole_l5_cross2")
l15c_1 <- load_slice("data/petiole2_A_sub2/",  "petiole_l15_cross1")
l15c_2 <- load_slice("data/petiole2_B_sub2/",  "petiole_l15_cross2")

# ----------------------------- #
# 2. Automatically split multi-section libraries #######
# ----------------------------- #

# Libraries the user identified as containing multiple separate cross sections:
# petiole_l1_cross, petiole_l2_cross, petiole_l4_cross1, petiole_l4_cross2, petiole_l5_cross2

# Tune eps per library if needed.
split_params <- tribble(
  ~base_name,            ~eps, ~minPts,
  "petiole_l1_cross",     500,   10,
  "petiole_l2_cross",     500,   10,
  "petiole_l4_cross1",    500,   10,
  "petiole_l4_cross2",    130,   8,
  "petiole_l5_cross2",    150,   10
)

raw_list <- list(
  petiole_l1_cross   = l1c,
  petiole_l2_cross   = l2c,
  petiole_l4_cross1  = l4c_1,
  petiole_l4_cross2  = l4c_2,
  petiole_l5_cross1  = l5c_1,
  petiole_l5_cross2  = l5c_2,
  petiole_l15_cross1 = l15c_1,
  petiole_l15_cross2 = l15c_2
)

# diagnostic plots before splitting
pdf("QC/dbscan_section_detection_checks.pdf", width = 12, height = 6)
for (i in seq_len(nrow(split_params))) {
  nm <- split_params$base_name[i]
  print(
    plot_dbscan_check(
      raw_list[[nm]],
      eps = split_params$eps[i],
      minPts = split_params$minPts[i],
      title = nm
    )
  )
}
dev.off()

split_out <- list()
split_tables <- list()

for (i in seq_len(nrow(split_params))) {
  nm <- split_params$base_name[i]
  tmp <- split_spatial_sections(
    seu = raw_list[[nm]],
    base_name = nm,
    eps = split_params$eps[i],
    minPts = split_params$minPts[i]
  )
  split_out[names(tmp)] <- tmp
  split_tables[[nm]] <- attr(tmp, "centroids")
}

split_table_all <- bind_rows(split_tables, .id = "parent_library")
write.csv(split_table_all, "petiole_split_sections_dbscan_section_centroids.csv")

# carry over the already single-section libraries unchanged
single_section_list <- list(
  petiole_l5_cross1 = l5c_1,
  petiole_l15_cross1 = l15c_1,
  petiole_l15_cross2 = l15c_2
)

for (nm in names(single_section_list)) {
  single_section_list[[nm]]$parent_library <- nm
  single_section_list[[nm]]$section_rank <- 1
  single_section_list[[nm]]$section_auto <- "section_1"
}

# merged list of all final sections
cross_list <- c(split_out, single_section_list)

# enforce a stable order
final_order <- names(cross_list)
write_lines(final_order, "tables/final_section_order.txt")

# ----------------------------- #
# 3. Merge objects ##############
# ----------------------------- #

sobj_cross <- merge(
  cross_list[[1]],
  y = cross_list[-1],
  add.cell.ids = names(cross_list),
  project = "Petiole_Cross_split"
)

qsave(sobj_cross, "saved_obj/sobj_cross_split_raw.qs")

# free memory
rm(l1c, l2c, l4c_1, l4c_2, l5c_1, l5c_2, l15c_1, l15c_2, raw_list, split_out, single_section_list)
gc()

# ----------------------------- #
# 4. QC + integration + clustering ########
# ----------------------------- #
### scaling + normalization + initial clustering 
sobj_cross <- process_spatial_integration(sobj_cross, resolution = 0.4, dims = 1:30)
qsave(sobj_cross, "saved_obj/sobj_cross_split_harmony_res0.4.qs")

# QC plots
plot1 <- VlnPlot_scCustom(
  sobj_cross,
  features = c("nCount_Spatial", "nFeature_Spatial"),
  group.by = "orig.ident",
  plot_median = TRUE
) + NoLegend()

plot2 <- SpatialFeaturePlot(
  sobj_cross,
  features = c("nCount_Spatial", "nFeature_Spatial"),
  crop = FALSE, ncol = 6
) + theme(legend.position = "right")

# ggsave("QC/QC_petiole_cross_split.2026.4.pdf", plot1 / plot2, width = 14, height = 14)
pdf("QC/QC_violin_petiole_cross_split.2026.4.pdf", width = 11, height = 7)
print(plot1)
dev.off()
pdf("QC/QC_spatial_petiole_cross_split.2026.4.pdf", width = 16, height = 14)
print(plot2)
dev.off()

sum_stats <- describeBy(sobj_cross@meta.data, group = sobj_cross@meta.data$orig.ident, mat = TRUE)
write.csv(sum_stats, "QC/QC_stats_beforeFilter_petiole_cross_split.csv", row.names = FALSE)

# Resolution sweep
resolutions <- c(0.3, 0.35, 0.4, 0.5, 0.6)
pdf("clustering/cross_resolution_optimization_split_L15.pdf", width = 15, height = 12)
for (res in resolutions) {
  message("Running clustering at resolution: ", res)
  sobj_tmp <- FindClusters(sobj_cross, resolution = res)
  p <- SpatialDimPlot(
    sobj_tmp,
    cells.highlight = CellsByIdentities(sobj_tmp),
    facet.highlight = TRUE,
    images = "petiole_l15_cross1",## petiole_l4_cross2_s3, petiole_l5_cross2_s1,petiole_l5_cross1
    ncol = 3,
    pt.size.factor = 4,
    alpha = 0.8
  ) +
    plot_annotation(
      title = paste0("Clustering resolution: ", res),
      theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    )
  print(p)
}
dev.off()


### Final resolution: 0.35 #######
sobj_cross <- FindClusters(sobj_cross, resolution = 0.35)

p1 <- DimPlot(sobj_cross, reduction = "umap", group.by = "orig.ident", alpha = 0.7) +
  ggtitle("By section")
p2 <- DimPlot(sobj_cross, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("By cluster")
p3 <- SpatialDimPlot(sobj_cross, group.by = "seurat_clusters", crop = FALSE, ncol = 4)
# final clustering report
pdf("clustering/petiole_cross_split_integration_res0.35.pdf", width = 14, height = 7)
p1+p2
p3
dev.off()

## Subcluster after 5. marker annotation (optional)  ####
## sub0 (phloem cambium)
sobj_cross_res0.35_sub0 <- FindSubCluster(sobj_cross, cluster = '0', "SCT_snn", subcluster.name = "sub0", resolution = 0.3, algorithm = 1)
Idents(sobj_cross_res0.35_sub0) <- sobj_cross_res0.35_sub0$sub0

cells_00 <- colnames(sobj_cross_res0.35_sub0)[sobj_cross_res0.35_sub0$sub0 == "0_0"]
cells_01 <- colnames(sobj_cross_res0.35_sub0)[sobj_cross_res0.35_sub0$sub0 == "0_1"]
SpatialDimPlot(
  sobj_cross_res0.35_sub0,
  cells.highlight = list(cells_00, cells_01),
  facet.highlight = TRUE,
  images = c("petiole_l15_cross1"),
  crop = TRUE)
DimPlot(sobj_cross_res0.35_sub0)
## sub5 (epidermal cortex)
sobj_cross_res0.35_sub05 <- FindSubCluster(sobj_cross_res0.35_sub0, cluster = '5', "SCT_snn", subcluster.name = "sub5", resolution = 0.1, algorithm = 1)
Idents(sobj_cross_res0.35_sub05) <- sobj_cross_res0.35_sub05$sub5

DimPlot(sobj_cross_res0.35_sub05)
cells_50 <- colnames(sobj_cross_res0.35_sub05)[sobj_cross_res0.35_sub05$sub5 == "5_0"]
cells_51 <- colnames(sobj_cross_res0.35_sub05)[sobj_cross_res0.35_sub05$sub5 == "5_1"]
# cells_52 <- colnames(sobj_cross_res0.35_sub05)[sobj_cross_res0.35_sub05$sub5 == "5_2"]
# cells_53 <- colnames(sobj_cross_res0.35_sub05)[sobj_cross_res0.35_sub05$sub5 == "5_3"]
SpatialDimPlot(
  sobj_cross_res0.35_sub05,
  cells.highlight = list(cells_50, cells_51),
  facet.highlight = TRUE,
  images = c("petiole_l5_cross1"),
  crop = TRUE)

### asign the final resolution
sobj_cross = sobj_cross_res0.35_sub05
sobj_cross$seurat_clusters <- Idents(sobj_cross) 

## Save the object with final_clustering -----------
qsave(sobj_cross, 'saved_obj/sobj_cross_harmony_res0.35_2026.4.qs')

sobj_cross <- qread("saved_obj/sobj_cross_harmony_res0.35_2026.4.qs")

## plot the proportion of cells in each cluster
plot_cluster_proportions(
  sobj = sobj_cross,
  order_vec = names(sobj_cross@images),
  out_pdf = "clustering/cross_integration_res0.35_ratio_labeled_ordered_split_subcluster.pdf"
)

## plot the split clusters for each section image:
cluster_cells <- CellsByIdentities(sobj_cross)
img_names <- names(sobj_cross@images)

pdf("clustering/clustering_split_cluster_cross_integration_res0.35_subcluster.pdf",
    width = 14, height = 15)

for (img in img_names) {
  p <- SpatialDimPlot(
    sobj_cross,
    cells.highlight = cluster_cells,
    images = img,
    facet.highlight = TRUE,
    alpha = 0.7,
    ncol = 4,
    pt.size.factor = 2,
    crop = FALSE,
    stroke = 0
  ) +
    plot_annotation(title = paste("Individual Cluster Spatial Distribution:", img))
  
  print(p)
}
dev.off()

# ----------------------------- #
# 5. Marker-based annotation ####
# ----------------------------- #
### Load the known marker list
mkr_list <- read.csv("updated_markerlist_poplar_2.8.26.csv")

# Keep only markers that exist
mkr_list_use <- mkr_list %>%
  filter(GeneID.v5 %in% rownames(sobj_cross))

# Heatmap for cluster interpretation
make_marker_heatmap(
  sobj = sobj_cross_res0.35_sub05,
  marker_df = mkr_list_use,
  out_pdf = "marker/marker_heatmap_split_celltype_petiole_cross_splitSections_subcluster.pdf"
)
#sobj_cross_res0.35_sub05

# provisional annotation by module scores
ann_res <- score_and_assign_celltypes(sobj_cross, mkr_list_use, assay = "SCT")
sobj_cross <- ann_res$object

write.csv(ann_res$cluster_scores, "tables/petiole_cross_cluster_marker_module_scores.csv")
write.csv(ann_res$top_assign, "tables/petiole_cross_cluster_top_predicted_celltypes.csv")
write.csv(ann_res$score_long, "tables/petiole_cross_cluster_marker_module_scores_long.csv")

# plots
plot_annotation_qc(
  sobj_cross,
  out_pdf = "clustering/petiole_cross_split_predicted_annotations_v2.pdf"
)

LinkedDimPlot(sobj_cross, image = 'petiole_l5_cross2_s1')

# Known-marker validation plots by broad tissue type
marker_groups <- list(
  Cortex = c("Cortex", "Ground meristem", "Endodermis"),
  Epidermis = c("Epidermis", "Guard cell", "Trichomes"),
  Phloem = c("Phloem", "Phloem Mother Cell", "Companion Cell", "Sieve Element",
             "Sieve element - companion cell complex"),
  Xylem = c("Xylem", "Vessel elements", "Ray cells"),
  Vascular = c("Cambium", "Vascular"),
  Other = c("Pith", "Cork", "Meristematic", "Shoot meristematic", "Proliferating")
)

img_groups <- split(names(sobj_cross@images), ceiling(seq_along(names(sobj_cross@images)) / 6))

for (grp in names(marker_groups)) {
  genes <- mkr_list_use %>%
    filter(celltype %in% marker_groups[[grp]]) %>%
    pull(GeneID.v5) %>%
    unique()
  
  if (length(genes) == 0) next
  
  gene_groups <- split(genes, ceiling(seq_along(genes) / 4))
  
  pdf(paste0("marker/known_marker_validation_", grp, "_splitSections.pdf"),
      width = 14, height = 18)
  
  for (gset in gene_groups) {
    for (imgs in img_groups) {
      print(
        SpatialFeaturePlot(
          sobj_cross,
          features = gset,
          images = imgs,
          ncol = 4,
          crop = FALSE,
          alpha = c(0.1, 1),
          pt.size.factor = 2
        ) +
          plot_annotation(title = paste(grp, "|", paste(gset, collapse = ", ")))
      )
    }
  }
  
  dev.off()
}

# optional: create a template file for manual final annotation after visual review
manual_map_template <- ann_res$top_assign %>%
  transmute(
    cluster_label,
    celltype_predicted = predicted_celltype,
    predicted_score,
    celltype_final = predicted_celltype
  ) %>%
  arrange(as.numeric(cluster_label))

write_csv(manual_map_template, "tables/manual_cluster_annotation_template_petiole_cross.csv")

# Plot the split clusters for sobj_cross using one of the images
pdf("clustering/petiole_cross_splitClusters_example.pdf", width = 12, height = 6)
print(
  SpatialDimPlot(
    sobj_cross,
    group.by = "seurat_clusters",
    images = names(sobj_cross@images)[1],
    label = FALSE,
    label.size = 4,
    pt.size.factor = 4
  )
)
dev.off()


# # If you edit tables/manual_cluster_annotation_template_petiole_cross.csv, rerun the block below:
# manual_map <- read_csv("tables/manual_cluster_annotation_template.csv", show_col_types = FALSE)
# sobj_cross <- rename_clusters_from_table(sobj_cross, manual_map)
# qsave(sobj_cross, "saved_obj/sobj_cross_split_harmony_res0.4_annotated.qs")

# For now, store predicted labels as celltypes_provisional
sobj_cross$celltypes_provisional <- sobj_cross$predicted_celltype
Idents(sobj_cross) <- sobj_cross$celltypes_provisional
SpatialDimPlot(sobj_cross, label = TRUE, images = "petiole_l5_cross2_s1")

pdf("clustering/petiole_cross_split_provisional_celltypes.pdf", width = 16, height = 8)
print(DimPlot(sobj_cross, reduction = "umap", group.by = "celltypes_provisional", label = TRUE))
print(SpatialDimPlot(sobj_cross, group.by = "celltypes_provisional", crop = FALSE, ncol = 4))
dev.off()

qsave(sobj_cross, "saved_obj/sobj_cross_split_harmony_res0.4_provisionalAnno.qs")


# ----------------------------- #
# 6. Rename cluster ##############
# ----------------------------- #
sobj_cross_ann <- RenameIdents(sobj_cross, '0_0' = 'Inner cortex','0_1'='Cortex','1'='Cortex', '2'='Epidermis',
                               '3'='Cortex','4'='Vasculature (Phloem)', '5_0'='Epidermis','5_1'='Epidermis', '6'='Pith (Xylem)', '7'='Cortex')
sobj_cross_ann$celltypes <- Idents(sobj_cross_ann)


## plot the proportion of cells in each cell type
plot_cluster_proportions(
  sobj = sobj_cross_ann,
  order_vec = names(sobj_cross_ann@images),
  out_pdf = "clustering/cross_integration_res0.35_ratio_labeled_ordered_split_subcluster_annotated.pdf"
)


## Plot UMAP and Spatial map of cell types
p1 <- DimPlot(sobj_cross_ann, reduction = "umap", group.by = "orig.ident", alpha = 0.7) 
p2 <- DimPlot(sobj_cross_ann, reduction = "umap", group.by = "celltypes", label = TRUE) 
  ggtitle("By cluster")
p3 <- SpatialDimPlot(sobj_cross_ann, group.by = "celltypes", crop = FALSE, ncol = 4)
p4 <- SpatialDimPlot(sobj_cross_ann, group.by = "celltypes", crop = TRUE, ncol = 4, pt.size.factor = 5, alpha = 0.9)
# final clustering report
pdf("clustering/petiole_cross_split_integration_res0.35_annotated.pdf", width = 15, height = 7)
p1+p2
p3
p4
dev.off()

### Save the annotated object #####
qsave(sobj_cross_ann, 'saved_obj/sobj_sp_petiole_cross_res0.35_anno_2026.4.qs')

sobj_cross_ann <- qread('saved_obj/sobj_sp_petiole_cross_res0.35_anno_2026.4.qs')
# ----------------------------- #

sobj_cross_ann <- atlas_annotations(sobj_cross_ann, "petiole_cross")
qs::qsave(sobj_cross_ann, "saved_obj/sobj_petiole_cross_res0.35_final_v1_2026.7.24.qs")

  },
  markers = function() {
sobj_cross_ann <- qs::qread(atlas_input("saved_obj/sobj_petiole_cross_res0.35_final_v1_2026.7.24.qs"))
Seurat::Idents(sobj_cross_ann) <- sobj_cross_ann$celltypes
# 7. Marker discovery.    ########
# ----------------------------- #

# Idents(sobj_cross_ann) <- sobj_cross_ann$celltypes
all_de_markers_petiole_cross <- FindAllMarkers(
  sobj_cross_ann,
  recorrect_umi = FALSE,
  test.use = "wilcox",
  logfc.threshold = 1,
  only.pos = TRUE,
  min.diff.pct = 0.2,
  min.pct = 0.1
)

write.csv(
  all_de_markers_petiole_cross,
  # "marker/sp_all_denovo_markers_petiole_cross_splitSections_res0.4_logfc1_20pctDif.csv",
  "marker/PETIOLE_CROSS_all_markers_by_celltype_PtXaOnly.csv", 
  row.names = FALSE
)

# Spatially variable genes
available_genes <- rownames(GetAssayData(sobj_cross, assay = "SCT", layer = "scale.data"))
top_features <- VariableFeatures(sobj_cross)[1:1000]
valid_features <- intersect(available_genes, top_features)

sobj_cross <- FindSpatiallyVariableFeatures(
  sobj_cross,
  assay = "SCT",
  slot = "scale.data",
  features = valid_features,
  selection.method = "moransi"
)

top_sv_features <- SpatiallyVariableFeatures(sobj_cross, method = "moransi")[1:15]

pdf("marker/spatially_variable_genes_moranI_petiole_cross_splitSections.pdf", width = 14, height = 12)
print(
  SpatialFeaturePlot(
    sobj_cross,
    features = top_sv_features,
    ncol = 4,
    alpha = c(0.1, 1),
    crop = FALSE
  )
)
dev.off()

qsave(sobj_cross, "saved_obj/sobj_cross_split_harmony_res0.4_provisionalAnno_SVG.qs")

message("Done. Main outputs:")
message(" - saved_obj/sobj_cross_split_harmony_res0.4_provisionalAnno.qs")
message(" - tables/manual_cluster_annotation_template_petiole_cross.csv")
message(" - clustering/petiole_cross_split_predicted_annotations.pdf")
message(" - marker/marker_heatmap_split_celltype_petiole_cross_splitSections.pdf")


# --------------------------------------- #

  },
  polarity = function() {
sobj_cross_ann <- qs::qread(atlas_input("saved_obj/sobj_petiole_cross_res0.35_final_v1_2026.7.24.qs"))
Seurat::Idents(sobj_cross_ann) <- sobj_cross_ann$celltypes
# 8. Cell type specific stage DEG ########
# --------------------------------------- #

DefaultAssay(sobj_cross_ann) <- "SCT"

stage_order <- c("L1", "L2", "L4", "L5", "L15")

sobj_cross_ann$stage <- case_when(
  str_detect(sobj_cross_ann$parent_library, "l15") ~ "L15",
  str_detect(sobj_cross_ann$parent_library, "l5")  ~ "L5",
  str_detect(sobj_cross_ann$parent_library, "l4")  ~ "L4",
  str_detect(sobj_cross_ann$parent_library, "l2")  ~ "L2",
  str_detect(sobj_cross_ann$parent_library, "l1")  ~ "L1",
  TRUE ~ NA_character_
)

sobj_cross_ann$stage <- factor(sobj_cross_ann$stage, levels = stage_order)

table(sobj_cross_ann$celltypes, sobj_cross_ann$stage)
table(sobj_cross_ann$parent_library, sobj_cross_ann$stage)

celltypes_to_test <- sort(unique(na.omit(sobj_cross_ann$celltypes)))

## 8.1 Pairwise stage DEGs within each cell type ----
run_stage_de <- function(obj, celltype_name, stage1, stage2, min_cells = 30) {
  sub <- subset(
    obj,
    subset = celltypes == celltype_name & stage %in% c(stage1, stage2)
  )
  
  if (ncol(sub) == 0) return(NULL)
  
  tab <- table(sub$stage)
  if (any(tab[c(stage1, stage2)] < min_cells)) return(NULL)
  
  Idents(sub) <- "stage"
  
  FindMarkers(
    sub,
    ident.1 = stage2,
    ident.2 = stage1,
    assay = "SCT",
    test.use = "wilcox",
    recorrect_umi = FALSE,
    logfc.threshold = 0.25,
    min.pct = 0.1
  ) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(
      celltype = celltype_name,
      stage1 = stage1,
      stage2 = stage2,
      contrast = paste0(stage2, "_vs_", stage1)
    )
}

stage_pairs <- list(
  c("L2", "L1"),
  c("L4", "L2"),
  c("L5", "L1"),
  c("L5", "L2"),
  c("L5", "L4"),
  c("L15", "L5"),
  c("L15", "L1")
)

de_stage_pairs_all <- purrr::map_dfr(celltypes_to_test, function(ct) {
  purrr::map_dfr(stage_pairs, function(sp) {
    run_stage_de(
      sobj_cross_ann,
      celltype_name = ct,
      stage1 = sp[2],
      stage2 = sp[1]
    )
  })
})

write.csv(
  de_stage_pairs_all,
  "tables/DE_petiole_cross_stage_by_celltype_all_pairs.csv",
  row.names = FALSE
)


## 8.2 Stage-specific DEGs: one stage vs all other stages within each cell type ----
run_stage_one_vs_all_de <- function(obj, celltype_name, stage_test, min_cells = 30) {
  other_stages <- setdiff(stage_order, stage_test)
  
  sub <- subset(
    obj,
    subset = celltypes == celltype_name & stage %in% c(stage_test, other_stages)
  )
  
  if (ncol(sub) == 0) return(NULL)
  
  tab <- table(sub$stage)
  n_test <- tab[stage_test]
  n_other <- sum(tab[other_stages])
  
  if (is.na(n_test) || n_test < min_cells || n_other < min_cells) return(NULL)
  
  Idents(sub) <- "stage"
  
  FindMarkers(
    sub,
    ident.1 = stage_test,
    ident.2 = other_stages,
    assay = "SCT",
    test.use = "wilcox",
    recorrect_umi = FALSE,
    logfc.threshold = 0.25,
    min.pct = 0.1
  ) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(
      celltype = celltype_name,
      stage = stage_test,
      contrast = paste0(stage_test, "_vs_all_other_stages"),
      direction = ifelse(avg_log2FC > 0, "stage_up", "other_stages_up")
    )
}

de_stage_one_vs_all <- purrr::map_dfr(celltypes_to_test, function(ct) {
  purrr::map_dfr(stage_order, function(st) {
    run_stage_one_vs_all_de(
      sobj_cross_ann,
      celltype_name = ct,
      stage_test = st
    )
  })
})

write.csv(
  de_stage_one_vs_all,
  "tables/DE_petiole_cross_stage_by_celltype_one_vs_all.csv",
  row.names = FALSE
)


## 8.3 Significant stage-specific genes for downstream GO ----
de_stage_one_vs_all_sig <- de_stage_one_vs_all %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  arrange(celltype, stage, p_val_adj)

write.csv(
  de_stage_one_vs_all_sig,
  "tables/DE_petiole_cross_stage_by_celltype_one_vs_all_sig.csv",
  row.names = FALSE
)
#
# de_stage_one_vs_all <- read.csv('tables/DE_petiole_cross_stage_by_celltype_one_vs_all.csv')
# de_stage_one_vs_all_sig <- read.csv('tables/DE_petiole_cross_stage_by_celltype_one_vs_all_sig.csv')

## 8.4 Summarize stage-specific DEGs ----
stage_order <- c("L1", "L2", "L4", "L5", "L15")

de_stage_one_vs_all_sig <- de_stage_one_vs_all %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
  mutate(stage = factor(trimws(as.character(stage)), levels = stage_order, ordered = TRUE)) %>%
  arrange(celltype, stage, p_val_adj)

write.csv(de_stage_one_vs_all_sig, "tables/DE_petiole_cross_stage_by_celltype_one_vs_all_sig.csv", row.names = FALSE)

stage_deg_summary <- de_stage_one_vs_all_sig %>%
  group_by(stage) %>%
  summarise(n_deg = n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(stage = factor(stage_order, levels = stage_order, ordered = TRUE), fill = list(n_deg = 0))

write.csv(stage_deg_summary, "tables/DE_petiole_cross_stage_specific_DEG_summary_by_stage.csv", row.names = FALSE)

p_stage_deg <- ggplot(stage_deg_summary, aes(x = n_deg, y = stage, group = 1)) +
  geom_path(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_discrete(limits = rev(stage_order)) +
  theme_classic() +
  labs(x = "Number of stage-specific DEGs", y = "Petiole stage")

ggsave("DE_petiole_cross_stage_specific_DEG_summary_by_stage_line.pdf", p_stage_deg, width = 5, height = 6)


p_stage_deg <- ggplot(stage_deg_summary, aes(x = n_deg, y = stage)) +
  geom_col(width = 0.7) +
  scale_y_discrete(limits = rev(stage_order)) +
  theme_classic() +
  labs(x = "Number of stage-specific DEGs", y = "Petiole stage")

ggsave("DE_petiole_cross_stage_specific_DEG_summary_by_stage.pdf", p_stage_deg, width = 5, height = 4)

## 8.5 Stage-specific DEG summary by cell type ----

stage_deg_by_celltype <- de_stage_one_vs_all_sig %>%
  group_by(celltype, stage) %>%
  summarise(n_deg = n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(celltype, stage = factor(stage_order, levels = stage_order, ordered = TRUE), fill = list(n_deg = 0)) %>%
  mutate(stage = factor(as.character(stage), levels = stage_order, ordered = TRUE))

write.csv(stage_deg_by_celltype, "tables/DE_petiole_cross_stage_specific_DEG_summary_by_celltype_stage.csv", row.names = FALSE)

p_stage_celltype_line <- ggplot(stage_deg_by_celltype, aes(x = n_deg, y = stage, group = celltype, color = celltype)) +
  geom_path(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_discrete(limits = rev(stage_order)) +
  theme_classic() +
  labs(x = "Number of stage-specific DEGs", y = "Petiole stage", color = "Cell type")

ggsave("DE_petiole_cross_stage_specific_DEG_summary_by_celltype_stage_line.pdf", p_stage_celltype_line, width = 6, height = 7)

## Heatmap 
p_stage_heatmap <- ggplot(stage_deg_by_celltype, aes(x = celltype, y = stage, fill = n_deg)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(colors = c("white", "#bdd7e7", "#6baed6", "#2171b5", "#08306b")) +
  scale_y_discrete(limits = rev(stage_order)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell type", y = "Petiole stage", fill = "DEGs")

ggsave("DE_petiole_cross_stage_specific_DEG_summary_by_celltype_stage_heatmap.pdf", p_stage_heatmap, width = 4, height = 4.5)
## --------------------------------------- #
# 9. Adaxial / abaxial DEG ###############
## --------------------------------------- #

## 9.1 helper functions ------------------

get_all_section_coords <- function(obj) {
  bind_rows(lapply(names(obj@images), function(img) {
    coords <- GetTissueCoordinates(obj[[img]])
    if (!"cell" %in% colnames(coords)) coords <- tibble::rownames_to_column(coords, "cell")
    coords$orig.ident <- img
    coords
  }))
}

get_section_meta_coords <- function(obj, img) {
  coords <- GetTissueCoordinates(obj[[img]])
  if (!"cell" %in% colnames(coords)) coords <- tibble::rownames_to_column(coords, "cell")
  
  obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    filter(orig.ident == img) %>%
    left_join(coords, by = "cell")
}

pick_line_points <- function(obj, img) {
  meta1 <- get_section_meta_coords(obj, img)
  
  plot(meta1$x, -meta1$y, asp = 1, pch = 16, cex = 0.8,
       xlab = "x", ylab = "-y", main = img)
  
  pts <- locator(2)
  
  tibble(
    orig.ident = img,
    x1 = pts$x[1], y1 = -pts$y[1],
    x2 = pts$x[2], y2 = -pts$y[2],
    adaxial_sign = 1
  )
}

assign_adaxial_abaxial_from_line <- function(obj, line_map) {
  coords <- get_all_section_coords(obj)
  
  meta2 <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    left_join(coords, by = c("cell", "orig.ident")) %>%
    left_join(line_map, by = "orig.ident") %>%
    mutate(
      side_val = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1),
      side_score = side_val * adaxial_sign,
      adaxial_abaxial = case_when(
        is.na(side_score) ~ NA_character_,
        side_score > 0 ~ "adaxial",
        TRUE ~ "abaxial"
      )
    )
  
  obj$side_score <- meta2$side_score[match(colnames(obj), meta2$cell)]
  obj$adaxial_abaxial <- meta2$adaxial_abaxial[match(colnames(obj), meta2$cell)]
  obj
}

assign_adaxial_abaxial_quarters <- function(obj, q = 0.25) {
  meta_q <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    group_by(orig.ident) %>%
    mutate(
      q_low = quantile(side_score, probs = q, na.rm = TRUE),
      q_high = quantile(side_score, probs = 1 - q, na.rm = TRUE),
      adaxial_abaxial_quarter = case_when(
        is.na(side_score) ~ NA_character_,
        side_score >= q_high ~ "adaxial_quarter",
        side_score <= q_low ~ "abaxial_quarter",
        TRUE ~ "middle"
      )
    ) %>%
    ungroup()
  
  obj$adaxial_abaxial_quarter <- meta_q$adaxial_abaxial_quarter[match(colnames(obj), meta_q$cell)]
  obj
}

run_side_de <- function(obj, ct, min_cells = 20) {
  sub <- subset(obj, subset = celltypes == ct & !is.na(adaxial_abaxial))
  tab <- table(sub$adaxial_abaxial)
  
  if (ncol(sub) == 0 || !all(c("adaxial", "abaxial") %in% names(tab)) ||
      any(tab[c("adaxial", "abaxial")] < min_cells)) return(NULL)
  
  Idents(sub) <- "adaxial_abaxial"
  
  FindMarkers(sub, ident.1 = "adaxial", ident.2 = "abaxial",
              assay = "SCT", test.use = "wilcox", recorrect_umi = FALSE,
              logfc.threshold = 0.25, min.pct = 0.1) %>%
    rownames_to_column("gene") %>%
    mutate(celltype = ct, comparison = "adaxial_vs_abaxial_all_stages")
}

run_side_stage_de <- function(obj, ct, st, min_cells = 15) {
  sub <- subset(obj, subset = celltypes == ct & stage == st & !is.na(adaxial_abaxial))
  tab <- table(sub$adaxial_abaxial)
  
  if (ncol(sub) == 0 || !all(c("adaxial", "abaxial") %in% names(tab)) ||
      any(tab[c("adaxial", "abaxial")] < min_cells)) return(NULL)
  
  Idents(sub) <- "adaxial_abaxial"
  
  FindMarkers(sub, ident.1 = "adaxial", ident.2 = "abaxial",
              assay = "SCT", test.use = "wilcox", recorrect_umi = FALSE,
              logfc.threshold = 0.25, min.pct = 0.1) %>%
    rownames_to_column("gene") %>%
    mutate(celltype = ct, stage = st, comparison = "adaxial_vs_abaxial")
}

run_side_quarter_de <- function(obj, ct, min_cells = 15) {
  sub <- subset(obj, subset = celltypes == ct & adaxial_abaxial_quarter %in% c("adaxial_quarter", "abaxial_quarter"))
  tab <- table(sub$adaxial_abaxial_quarter)
  
  if (ncol(sub) == 0 || !all(c("adaxial_quarter", "abaxial_quarter") %in% names(tab)) ||
      any(tab[c("adaxial_quarter", "abaxial_quarter")] < min_cells)) return(NULL)
  
  Idents(sub) <- "adaxial_abaxial_quarter"
  
  FindMarkers(sub, ident.1 = "adaxial_quarter", ident.2 = "abaxial_quarter",
              assay = "SCT", test.use = "wilcox", recorrect_umi = FALSE,
              logfc.threshold = 0.25, min.pct = 0.1) %>%
    rownames_to_column("gene") %>%
    mutate(celltype = ct, comparison = "adaxial_quarter_vs_abaxial_quarter_all_stages")
}

run_side_stage_quarter_de <- function(obj, ct, st, min_cells = 10) {
  sub <- subset(obj, subset = celltypes == ct & stage == st &
                  adaxial_abaxial_quarter %in% c("adaxial_quarter", "abaxial_quarter"))
  tab <- table(sub$adaxial_abaxial_quarter)
  
  if (ncol(sub) == 0 || !all(c("adaxial_quarter", "abaxial_quarter") %in% names(tab)) ||
      any(tab[c("adaxial_quarter", "abaxial_quarter")] < min_cells)) return(NULL)
  
  Idents(sub) <- "adaxial_abaxial_quarter"
  
  FindMarkers(sub, ident.1 = "adaxial_quarter", ident.2 = "abaxial_quarter",
              assay = "SCT", test.use = "wilcox", recorrect_umi = FALSE,
              logfc.threshold = 0.25, min.pct = 0.1) %>%
    rownames_to_column("gene") %>%
    mutate(celltype = ct, stage = st, comparison = "adaxial_quarter_vs_abaxial_quarter")
}

summarize_side_de <- function(de_df) {
  de_df %>%
    mutate(sig = p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
    group_by(celltype, stage) %>%
    summarise(
      n_deg = sum(sig, na.rm = TRUE),
      n_adaxial_up = sum(sig & avg_log2FC > 0, na.rm = TRUE),
      n_abaxial_up = sum(sig & avg_log2FC < 0, na.rm = TRUE),
      .groups = "drop"
    )
}

## 9.2 manually draw boundary line for each section ----

img_names <- names(sobj_cross_ann@images)

line_map <- bind_rows(lapply(img_names, function(img) {
  message("Click 2 points for: ", img)
  pick_line_points(sobj_cross_ann, img)
}))

write.csv(line_map, "tables/adaxial_abaxial_line_map_raw.csv", row.names = FALSE)

line_map <- read.csv("tables/adaxial_abaxial_line_map_raw.csv")

## optional: flip any wrongly assigned section and rerun assignment later
# line_map$adaxial_sign[line_map$orig.ident == "petiole_l4_cross2_s1"] <- -1

## 9.3 assign adaxial / abaxial labels ----

sobj_cross_ann <- assign_adaxial_abaxial_from_line(sobj_cross_ann, line_map)
sobj_cross_ann <- assign_adaxial_abaxial_quarters(sobj_cross_ann, q = 0.3)

img_names <- names(sobj_cross_ann@images)
pdf("clustering/adaxial_abaxial_check_all_sections.pdf", width = 8, height = 8)
for (img in img_names) {
  print(SpatialDimPlot(sobj_cross_ann, group.by = "adaxial_abaxial",
                       images = img, crop = FALSE) + ggtitle(img))
}
dev.off()

pdf("clustering/adaxial_abaxial_quarter_check_all_sections.pdf", width = 8, height = 8)
for (img in img_names) {
  print(SpatialDimPlot(sobj_cross_ann, group.by = "adaxial_abaxial_quarter",
                       images = img, crop = FALSE) + ggtitle(img))
}
dev.off()

write.csv(line_map, "tables/adaxial_abaxial_line_map_final.csv", row.names = FALSE)

table(sobj_cross_ann$adaxial_abaxial)
table(sobj_cross_ann$adaxial_abaxial_quarter)

## 9.4 global adaxial vs abaxial DE by cell type ----

stage_order <- c("L1", "L2", "L4", "L5", "L15")
celltypes_use <- sort(unique(sobj_cross_ann$celltypes))
stages_use <- stage_order[stage_order %in% unique(as.character(sobj_cross_ann$stage))]

de_side_all <- map(celltypes_use, ~ run_side_de(sobj_cross_ann, .x)) %>%
  bind_rows()

write.csv(de_side_all, "tables/DE_adaxial_abaxial_by_celltype_all_stages.csv", row.names = FALSE)

## 9.5 stage-specific adaxial vs abaxial DE ----

de_side_stage_all <- crossing(celltype = celltypes_use, stage = stages_use) %>%
  mutate(res = map2(celltype, stage, ~ run_side_stage_de(sobj_cross_ann, .x, .y))) %>%
  pull(res) %>%
  bind_rows()

write.csv(de_side_stage_all, "tables/DE_adaxial_abaxial_by_celltype_each_stage.csv", row.names = FALSE)
# de_side_stage_all <- read.csv("tables/DE_adaxial_abaxial_by_celltype_each_stage.csv")

## 9.6 quarter-based adaxial vs abaxial DE ----

de_side_quarter_all <- map(celltypes_use, ~ run_side_quarter_de(sobj_cross_ann, .x)) %>%
  bind_rows()

write.csv(de_side_quarter_all, "tables/DE_adaxial_abaxial_quarter_by_celltype_all_stages.csv", row.names = FALSE)

de_side_stage_quarter_all <- crossing(celltype = celltypes_use, stage = stages_use) %>%
  mutate(res = map2(celltype, stage, ~ run_side_stage_quarter_de(sobj_cross_ann, .x, .y))) %>%
  pull(res) %>%
  bind_rows()

write.csv(de_side_stage_quarter_all, "tables/DE_adaxial_abaxial_quarter_by_celltype_each_stage.csv", row.names = FALSE)

## 9.7 summary: when does asymmetry start? ----

de_side_stage_summary <- summarize_side_de(de_side_stage_all)
de_side_stage_quarter_summary <- summarize_side_de(de_side_stage_quarter_all)

write.csv(de_side_stage_summary, "tables/DE_adaxial_abaxial_by_celltype_stage_summary.csv", row.names = FALSE)
write.csv(de_side_stage_quarter_summary, "tables/DE_adaxial_abaxial_quarter_by_celltype_stage_summary.csv", row.names = FALSE)

first_asymmetry <- de_side_stage_summary %>%
  mutate(stage = factor(stage, levels = stage_order, ordered = TRUE)) %>%
  filter(n_deg > 0) %>%
  arrange(celltype, stage) %>%
  group_by(celltype) %>%
  slice_head(n = 1) %>%
  ungroup()

first_asymmetry_quarter <- de_side_stage_quarter_summary %>%
  mutate(stage = factor(stage, levels = stage_order, ordered = TRUE)) %>%
  filter(n_deg > 0) %>%
  arrange(celltype, stage) %>%
  group_by(celltype) %>%
  slice_head(n = 1) %>%
  ungroup()

write.csv(first_asymmetry, "tables/DE_adaxial_abaxial_first_stage_per_celltype.csv", row.names = FALSE)
write.csv(first_asymmetry_quarter, "tables/DE_adaxial_abaxial_quarter_first_stage_per_celltype.csv", row.names = FALSE)

## 9.8 plot DEG counts by stage ----

de_side_stage_summary <- de_side_stage_summary %>%
  mutate(stage = factor(trimws(as.character(stage)), levels = stage_order, ordered = TRUE))

de_side_stage_quarter_summary <- de_side_stage_quarter_summary %>%
  mutate(stage = factor(trimws(as.character(stage)), levels = stage_order, ordered = TRUE))

p_deg <- ggplot(de_side_stage_summary, aes(stage, n_deg, group = celltype, color = celltype)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_discrete(limits = stage_order, drop = FALSE) +
  theme_classic() +
  labs(x = "Stage", y = "Number of adaxial/abaxial DE genes", color = "Cell type")

ggsave("DE_adaxial_abaxial_by_celltype_stage_summary.pdf", p_deg, width = 12, height = 4)

p_deg_quarter <- ggplot(de_side_stage_quarter_summary, aes(stage, n_deg, group = celltype, color = celltype)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_discrete(limits = stage_order, drop = FALSE) +
  theme_classic() +
  labs(x = "Stage", y = "Number of upper/lower quarter DE genes", color = "Cell type")

ggsave("DE_adaxial_abaxial_quarter_by_celltype_stage_summary.pdf", p_deg_quarter, width = 12, height = 4)

de_compare_summary <- bind_rows(
  de_side_stage_summary %>% mutate(method = "split_half"),
  de_side_stage_quarter_summary %>% mutate(method = "upper_lower_quarter")
)

write.csv(de_compare_summary, "tables/DE_adaxial_abaxial_half_vs_quarter_summary.csv", row.names = FALSE)

p_compare_deg <- ggplot(de_compare_summary, aes(stage, n_deg, group = interaction(celltype, method),
                                                color = celltype, linetype = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_discrete(limits = stage_order, drop = FALSE) +
  theme_classic() +
  labs(x = "Stage", y = "Number of DE genes", color = "Cell type", linetype = "Comparison")

ggsave("DE_adaxial_abaxial_half_vs_quarter_summary.pdf", p_compare_deg, width = 13, height = 5)

### 9.8b vertical version, original half-line split ----

de_side_stage_summary_v2 <- de_side_stage_summary %>%
  arrange(celltype, stage)

p_deg_v2 <- ggplot(de_side_stage_summary_v2, aes(x = n_deg, y = stage, group = celltype, color = celltype)) +
  geom_path(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_discrete(limits = rev(stage_order)) +
  theme_classic() +
  labs(x = "Number of adaxial/abaxial DE genes", y = "Stage", color = "Cell type")

ggsave("DE_adaxial_abaxial_by_celltype_stage_summary_v2_vertical.pdf", p_deg_v2, width = 5, height = 8)

### 9.8c vertical version, quarter split ----

de_side_stage_quarter_summary_v2 <- de_side_stage_quarter_summary %>%
  arrange(celltype, stage)

p_deg_quarter_v2 <- ggplot(de_side_stage_quarter_summary_v2, aes(x = n_deg, y = stage, group = celltype, color = celltype)) +
  geom_path(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_y_discrete(limits = rev(stage_order)) +
  theme_classic() +
  labs(x = "Number of upper/lower quarter DE genes", y = "Stage", color = "Cell type")

ggsave("DE_adaxial_abaxial_quarter_by_celltype_stage_summary_v2_vertical.pdf", p_deg_quarter_v2, width = 5, height = 8)

## 9.9 export significant candidates ----

candidate_deg_side <- de_side_stage_all %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)

candidate_deg_side_quarter <- de_side_stage_quarter_all %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25)

write.csv(candidate_deg_side, "DE_adaxial_abaxial_significant_candidates.csv", row.names = FALSE)
write.csv(candidate_deg_side_quarter, "DE_adaxial_abaxial_quarter_significant_candidates.csv", row.names = FALSE)

## 9.10 compare candidate overlap between half-line and quarter methods ----

candidate_overlap <- inner_join(
  candidate_deg_side %>% select(gene, celltype, stage) %>% distinct() %>% mutate(half_line = TRUE),
  candidate_deg_side_quarter %>% select(gene, celltype, stage) %>% distinct() %>% mutate(quarter = TRUE),
  by = c("gene", "celltype", "stage")
)

candidate_union <- full_join(
  candidate_deg_side %>% select(gene, celltype, stage) %>% distinct() %>% mutate(half_line = TRUE),
  candidate_deg_side_quarter %>% select(gene, celltype, stage) %>% distinct() %>% mutate(quarter = TRUE),
  by = c("gene", "celltype", "stage")
) %>%
  mutate(
    half_line = ifelse(is.na(half_line), FALSE, half_line),
    quarter = ifelse(is.na(quarter), FALSE, quarter)
  )

write.csv(candidate_overlap, "tables/DE_adaxial_abaxial_half_quarter_overlap_candidates.csv", row.names = FALSE)
write.csv(candidate_union, "tables/DE_adaxial_abaxial_half_quarter_union_candidates.csv", row.names = FALSE)

## 9.11 optional: plot selected genes across stages ----

genes_to_plot <- c("PtXaAlbH.01G087200.v5.1", "PtXaTreH.06G108400.v5.1", "PtXaAlbH.14G011500.v5.1")  # replace

pdf("marker/adaxial_abaxial_candidate_genes_by_stage.pdf", width = 14, height = 10)
for (st in stage_order) {
  sub <- subset(sobj_cross_ann, subset = stage == st)
  feats <- intersect(genes_to_plot, rownames(sub))
  if (length(feats) == 0) next
  
  print(SpatialFeaturePlot(sub, features = feats, crop = FALSE, ncol = 3) +
          ggtitle(paste("Stage:", st)))
}
dev.off()

## save the object with updated meta-data (stages, sides)
qsave(sobj_cross_ann, "saved_obj/sobj_petiole_cross_res0.35_final_v1_2026.7.24.qs")



## 9.12 Filter asymmetric DEGs for true side-enriched genes --------

# Goal:
# Remove genes that are likely only broad cell-type markers.
# Prioritize genes that show adaxial/abaxial spatial enrichment within the same cell type/stage.

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)

stage_order <- c("L1", "L2", "L4", "L5", "L15")


#### 9.12.1 load candidate DEGs from half and quarter methods #########

deg_half <- read.csv("DE_adaxial_abaxial_significant_candidates.csv") %>%
  mutate(method = "split_half")

deg_quarter <- read.csv("DE_adaxial_abaxial_quarter_significant_candidates.csv") %>%
  mutate(method = "upper_lower_quarter")

deg_side_candidates <- bind_rows(deg_half, deg_quarter) %>%
  mutate(stage = factor(trimws(as.character(stage)), levels = stage_order, ordered = TRUE),
         side_direction = ifelse(avg_log2FC > 0, "Adaxial up", "Abaxial up"))

write.csv(deg_side_candidates, "tables/DE_adaxial_abaxial_all_candidate_methods.csv", row.names = FALSE)


### 9.12.2 remove broad cell-type markers ####

marker_file <- "marker/PETIOLE_CROSS_all_markers_by_celltype_PtXaOnly.csv"

if (file.exists(marker_file)) {
  celltype_markers <- read.csv(marker_file) %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 100, with_ties = FALSE) %>%
    ungroup() %>%
    pull(gene) %>%
    unique()
} else {
  message("Celltype marker file not found. Recalculate markers from sobj_cross_ann.")
  
  Idents(sobj_cross_ann) <- sobj_cross_ann$celltypes
  celltype_marker_df <- FindAllMarkers(
    sobj_cross_ann,
    assay = "SCT",
    recorrect_umi = FALSE,
    test.use = "wilcox",
    only.pos = TRUE,
    logfc.threshold = 0.5,
    min.pct = 0.1
  )
  
  write.csv(celltype_marker_df, "marker/petiole_cross_all_markers_by_celltype.csv", row.names = FALSE)
  
  celltype_markers <- celltype_marker_df %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = 100, with_ties = FALSE) %>%
    ungroup() %>%
    pull(gene) %>%
    unique()
}

deg_side_candidates_filtered <- deg_side_candidates %>%
  mutate(is_celltype_marker = gene %in% celltype_markers) %>%
  filter(!is_celltype_marker)

write.csv(deg_side_candidates_filtered,
          "tables/DE_adaxial_abaxial_candidates_excluding_celltype_markers.csv",
          row.names = FALSE)

#### 9.12.3 test gene expression correlation with adaxial-abaxial axis ######

test_side_score_cor <- function(obj, gene, ct, st, min_cells = 20) {
  md <- obj@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(stage_chr = as.character(stage),
           celltypes_chr = as.character(celltypes))
  
  cells_use <- md %>%
    filter(celltypes_chr == as.character(ct),
           stage_chr == as.character(st),
           !is.na(side_score)) %>%
    pull(cell)
  
  if (length(cells_use) < min_cells || !gene %in% rownames(obj)) return(NULL)
  
  expr <- as.numeric(GetAssayData(obj, assay = "SCT", layer = "data")[gene, cells_use])
  side <- as.numeric(md$side_score[match(cells_use, md$cell)])
  
  if (sd(expr, na.rm = TRUE) == 0 || sd(side, na.rm = TRUE) == 0) return(NULL)
  
  ct_res <- suppressWarnings(cor.test(expr, side, method = "spearman"))
  
  tibble(gene = gene, celltype = as.character(ct), stage = as.character(st),
         rho_side_score = unname(ct_res$estimate),
         p_side_score = ct_res$p.value, n_cells = length(cells_use))
}

side_score_cor_res <- deg_side_candidates_filtered %>%
  distinct(gene, celltype, stage) %>%
  mutate(res = pmap(list(gene, celltype, stage),
                    ~ test_side_score_cor(sobj_cross_ann, ..1, ..2, ..3))) %>%
  pull(res) %>%
  bind_rows() %>%
  mutate(p_adj_side_score = p.adjust(p_side_score, method = "BH"))

write.csv(side_score_cor_res,
          "tables/DE_adaxial_abaxial_candidates_side_score_correlation.csv",
          row.names = FALSE)

## rho_side_score > 0  = higher toward adaxial side
## rho_side_score < 0  = higher toward abaxial side

deg_side_prioritized <- deg_side_candidates_filtered %>%
  left_join(side_score_cor_res,
            by = c("gene", "celltype", "stage")) %>%
  mutate(
    direction_consistent = case_when(
      side_direction == "Adaxial up" & rho_side_score > 0 ~ TRUE,
      side_direction == "Abaxial up" & rho_side_score < 0 ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(direction_consistent, p_adj_side_score < 0.05)

write.csv(deg_side_prioritized,
          "tables/DE_adaxial_abaxial_prioritized_spatially_enriched_genes.csv",
          row.names = FALSE)

#### 9.12.4 section-level reproducibility of side enrichment ######

test_gene_side_by_section <- function(obj, gene, ct, st, min_cells = 8) {
  md <- obj@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(stage_chr = as.character(stage),
           celltypes_chr = as.character(celltypes))
  
  cells_use <- md %>%
    filter(celltypes_chr == as.character(ct),
           stage_chr == as.character(st),
           !is.na(adaxial_abaxial)) %>%
    pull(cell)
  
  if (length(cells_use) == 0 || !gene %in% rownames(obj)) return(NULL)
  
  map_dfr(unique(md$orig.ident[match(cells_use, md$cell)]), function(sec) {
    cells_sec <- md %>%
      filter(cell %in% cells_use, orig.ident == sec) %>%
      pull(cell)
    
    md_sec <- md[match(cells_sec, md$cell), ]
    tab <- table(md_sec$adaxial_abaxial)
    
    if (!all(c("adaxial", "abaxial") %in% names(tab)) ||
        any(tab[c("adaxial", "abaxial")] < min_cells)) return(NULL)
    
    expr <- as.numeric(GetAssayData(obj, assay = "SCT", layer = "data")[gene, cells_sec])
    
    tibble(
      gene = gene,
      celltype = as.character(ct),
      stage = as.character(st),
      orig.ident = sec,
      mean_adaxial = mean(expr[md_sec$adaxial_abaxial == "adaxial"], na.rm = TRUE),
      mean_abaxial = mean(expr[md_sec$adaxial_abaxial == "abaxial"], na.rm = TRUE),
      diff_adaxial_abaxial = mean_adaxial - mean_abaxial
    )
  })
}

side_section_res <- deg_side_prioritized %>%
  distinct(gene, celltype, stage) %>%
  mutate(gene = as.character(gene),
         celltype = as.character(celltype),
         stage = as.character(stage)) %>%
  mutate(res = pmap(list(gene, celltype, stage),
                    ~ test_gene_side_by_section(sobj_cross_ann, ..1, ..2, ..3))) %>%
  pull(res) %>%
  bind_rows()

write.csv(side_section_res,
          "tables/DE_adaxial_abaxial_prioritized_section_reproducibility.csv",
          row.names = FALSE)

side_section_summary <- side_section_res %>%
  group_by(gene, celltype, stage) %>%
  summarise(
    n_sections = n(),
    n_sections_adaxial_up = sum(diff_adaxial_abaxial > 0, na.rm = TRUE),
    n_sections_abaxial_up = sum(diff_adaxial_abaxial < 0, na.rm = TRUE),
    mean_diff = mean(diff_adaxial_abaxial, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(side_section_summary,
          "tables/DE_adaxial_abaxial_prioritized_section_reproducibility_summary.csv",
          row.names = FALSE)

final_spatial_side_candidates <- deg_side_prioritized %>%
  left_join(side_section_summary, by = c("gene", "celltype", "stage")) %>%
  filter(
    n_sections >= 2,
    (side_direction == "Adaxial up" & n_sections_adaxial_up >= 2) |
      (side_direction == "Abaxial up" & n_sections_abaxial_up >= 2)
  )

write.csv(final_spatial_side_candidates,
          "tables/DE_adaxial_abaxial_final_spatial_side_candidates.csv",
          row.names = FALSE)

#### 9.12.5 plot final candidate genes within celltype/stage ####

## Select the top canditates that showed spatially asymmetric expression
n_top <- 50
top_plot_candidates <- final_spatial_side_candidates %>%
  arrange(p_adj_side_score, p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = n_top)

write.csv(top_plot_candidates, "tables/DE_adaxial_abaxial_top50_plot_candidates.csv", row.names = FALSE)

md <- sobj_cross_ann@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(stage_chr = as.character(stage), celltypes_chr = as.character(celltypes))


# plot one candidate manually: change i to any row number in top_plot_candidates
i <- 1
gene_i <- top_plot_candidates$gene[i]
ct_i <- top_plot_candidates$celltype[i]
st_i <- as.character(top_plot_candidates$stage[i])

cells_use <- md %>%
  filter(celltypes_chr == as.character(ct_i), stage_chr == as.character(st_i)) %>%
  pull(cell)

sub <- subset(sobj_cross_ann, cells = cells_use)
img_use <- names(sub@images)

p1 <- SpatialFeaturePlot(sub, features = gene_i, images = img_use, crop = FALSE, ncol = 2,
                         min.cutoff = "q05", max.cutoff = "q95") +
  ggtitle(paste(gene_i, ct_i, st_i, sep = " | "))

p2 <- VlnPlot(sub, features = gene_i, group.by = "adaxial_abaxial", pt.size = 0.1) +
  ggtitle(paste(gene_i, "Half split"))

p3 <- VlnPlot(sub, features = gene_i, group.by = "adaxial_abaxial_quarter", pt.size = 0.1) +
  ggtitle(paste(gene_i, "Quarter split"))

p1
p2
p3


##### plot selected genes of interest####

## SAUR-like auxin-responsive protein
g_for_plt <- c('PtXaTreH.06G108400.v5.1', 'PtXaTreH.09G107100.v5.1', 'PtXaTreH.01G139400.v5.1', 'PtXaTreH.05G167000.v5.1', 'PtXaAlbH.05G168000.v5.1')
## SAUR-like auxin-responsive protein: c('PtXaTreH.06G108400.v5.1', 'PtXaTreH.09G107100.v5.1', 'PtXaTreH.01G139400.v5.1')
## AUX/IAA transcriptional regulator family: c('PtXaTreH.05G167000.v5.1' , 'PtXaAlbH.05G168000.v5.1')
g_for_plt_df <- tibble(
  gene = g_for_plt,
  title = c("SAUR-like auxin-responsive protein", "SAUR-like auxin-responsive protein", "SAUR-like auxin-responsive protein", 'AUX/IAA family', 'AUX/IAA family'),
  celltype = c("Cortex", "Cortex", "Cortex","Cortex","Epidermis" ),  # edit
  stage = c("L4", "L5", "L5", "L5", "L1")                 # edit
)

g_for_plt <- c('PtXaAlbH.07G007900.v5.1', 'PtXaTreH.04G015300.v5.1', 'PtXaTreH.16G078300.v5.1', 'PtXaTreH.01G057400.v5.1', 'PtXaTreH.04G098000.v5.1', 'PtXaAlbH.04G097300.v5.1', 'PtXaTreH.04G098000.v5.1')
## cell-wall remodeling: c('PtXaAlbH.07G007900.v5.1', 'PtXaTreH.04G015300.v5.1', 'PtXaTreH.16G078300.v5.1', 'PtXaTreH.01G057400.v5.1', 'PtXaTreH.04G098000.v5.1', 'PtXaAlbH.04G097300.v5.1', 'PtXaTreH.04G098000.v5.1')
g_for_plt_df <- tibble(
  gene = g_for_plt,
  title = c("xyloglucan hydrolase", "xyloglucan hydrolase", "xyloglucan hydrolase", 'xyloglucan hydrolase', 'expansin A5','expansin A5','expansin A5'),
  celltype = c("Cortex", "Cortex", "Cortex","Cortex","Cortex", "Cortex", "Cortex"),  # edit
  stage = c("L4", "L4", "L4", "L4", "L4","L4","L5"))

md <- sobj_cross_ann@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(stage_chr = as.character(stage), celltypes_chr = as.character(celltypes))

pdf("marker/selected_cellwall_candidate_genes_spatial_pages.pdf", width = 12, height = 12)
for (i in seq_len(nrow(g_for_plt_df))) {
  gene_i <- g_for_plt_df$gene[i]; ct_i <- g_for_plt_df$celltype[i]; st_i <- as.character(g_for_plt_df$stage[i]); title_i <- g_for_plt_df$title[i]
  cells_use <- md %>% filter(celltypes_chr == as.character(ct_i), stage_chr == as.character(st_i)) %>% pull(cell)
  if (length(cells_use) == 0 || !gene_i %in% rownames(sobj_cross_ann)) next
  sub <- subset(sobj_cross_ann, cells = cells_use); img_use <- names(sub@images)
  p1 <- SpatialFeaturePlot(sub, features = gene_i, images = img_use, crop = FALSE, ncol = 2,
                           min.cutoff = "q05", max.cutoff = "q95") +
    ggtitle(paste(title_i, gene_i, ct_i, st_i, sep = " | "))
  print(p1)
}
dev.off()

pdf("marker/selected_auxin_candidate_genes_violin_pages.pdf", width = 4, height = 6)
for (i in seq_len(nrow(g_for_plt_df))) {
  gene_i <- g_for_plt_df$gene[i]; ct_i <- g_for_plt_df$celltype[i]; st_i <- as.character(g_for_plt_df$stage[i]); title_i <- g_for_plt_df$title[i]
  cells_use <- md %>% filter(celltypes_chr == as.character(ct_i), stage_chr == as.character(st_i)) %>% pull(cell)
  if (length(cells_use) == 0 || !gene_i %in% rownames(sobj_cross_ann)) next
  sub <- subset(sobj_cross_ann, cells = cells_use)
  p2 <- VlnPlot(sub, features = gene_i, group.by = "adaxial_abaxial", pt.size = 0.05) + theme_classic(base_size = 10) + theme(legend.position = "none")
  p3 <- VlnPlot(sub, features = gene_i, group.by = "adaxial_abaxial_quarter", pt.size = 0.05) + theme_classic(base_size = 10) + theme(legend.position = "none")
  print(p2 + p3)
}
dev.off()


## Plot all the top n candidate genes in top_plot_candidates
# spatial plots: larger page
pdf("marker/side_candidate_top_genes_spatial_pages.pdf", width = 12, height = 12)

for (i in seq_len(nrow(top_plot_candidates))) {
  gene_i <- top_plot_candidates$gene[i]; ct_i <- top_plot_candidates$celltype[i]; st_i <- as.character(top_plot_candidates$stage[i])
  cells_use <- md %>% filter(celltypes_chr == as.character(ct_i), stage_chr == as.character(st_i)) %>% pull(cell)
  if (length(cells_use) == 0 || !gene_i %in% rownames(sobj_cross_ann)) next
  sub <- subset(sobj_cross_ann, cells = cells_use); img_use <- names(sub@images)
  
  p1 <- SpatialFeaturePlot(sub, features = gene_i, images = img_use, crop = FALSE, ncol = 2,
                           min.cutoff = "q05", max.cutoff = "q95") +
  patchwork::plot_annotation(title = paste0("Rank ", i, " | ", gene_i, " | ", ct_i, " | ", st_i))
  print(p1)
}

dev.off()

# half split violin: compact page
# both violin plots on one compact page
pdf("marker/side_candidate_top_genes_violin_combined_pages.pdf", width = 4, height = 6)

for (i in seq_len(nrow(top_plot_candidates))) {
  gene_i <- top_plot_candidates$gene[i]; ct_i <- top_plot_candidates$celltype[i]; st_i <- as.character(top_plot_candidates$stage[i])
  cells_use <- md %>% filter(celltypes_chr == as.character(ct_i), stage_chr == as.character(st_i)) %>% pull(cell)
  if (length(cells_use) == 0 || !gene_i %in% rownames(sobj_cross_ann)) next
  sub <- subset(sobj_cross_ann, cells = cells_use)
  
  ttl <- paste0("Rank ", i, " | ", gene_i, " | ", ct_i, " | ", st_i)
  
  p2 <- VlnPlot(sub, features = gene_i, group.by = "adaxial_abaxial", pt.size = 0.05) +
    ggtitle(paste0(ttl, " | Half split")) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(size = 8), axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
  
  p3 <- VlnPlot(sub, features = gene_i, group.by = "adaxial_abaxial_quarter", pt.size = 0.05) +
    ggtitle(paste0(ttl, " | Quarter split")) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(size = 8), axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
  
  print(p2 + p3)
}

dev.off()

### 9.12.6 half-split DEG recurrence, overlap, and temporal pattern modules ####
## Purpose:
##   1) find genes exclusive to one celltype/stage/side contrast versus genes
##      recurrent across multiple petiole domains or stages;
##   2) provide an UpSet-style overlap view for half-split DEGs;
##   3) cluster recurrent side-biased genes by developmental expression pattern
##      to identify early, L4/L5-peak, late, and consistently high programs.

dir.create("overlap", showWarnings = FALSE, recursive = TRUE)
dir.create("mfuzz", showWarnings = FALSE, recursive = TRUE)

stage_order_overlap <- if (exists("stage_order")) stage_order else c("L1", "L2", "L4", "L5", "L15")
overlap_base_family <- "sans"

save_overlap_pdf <- function(plot, filename, width, height) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  tmp_file <- paste0(filename, ".tmp")
  if (file.exists(tmp_file)) unlink(tmp_file)
  dev_before <- grDevices::dev.cur()
  ok <- FALSE
  grDevices::cairo_pdf(tmp_file, width = width, height = height,
                       family = overlap_base_family, onefile = TRUE)
  tryCatch(
    {
      print(plot)
      ok <- TRUE
    },
    finally = {
      if (grDevices::dev.cur() != dev_before) grDevices::dev.off()
    }
  )
  if (ok) {
    if (file.exists(filename)) unlink(filename)
    file.rename(tmp_file, filename)
  }
  invisible(filename)
}

## 9.12.6a recurrence table from half-split significant DEGs ----
de_half_sig_for_overlap <- de_side_stage_all %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
  mutate(
    stage = factor(as.character(stage), levels = stage_order_overlap, ordered = TRUE),
    side_direction = ifelse(avg_log2FC > 0, "Adaxial up", "Abaxial up"),
    set_id = paste(celltype, stage, side_direction, sep = "__"),
    set_id_clean = make.names(set_id)
  )

write.csv(
  de_half_sig_for_overlap,
  "tables/DE_adaxial_abaxial_half_split_sig_for_overlap.csv",
  row.names = FALSE
)

half_split_set_key <- de_half_sig_for_overlap %>%
  distinct(set_id, set_id_clean, celltype, stage, side_direction) %>%
  arrange(celltype, stage, side_direction)

write.csv(
  half_split_set_key,
  "tables/DE_adaxial_abaxial_half_split_overlap_set_key.csv",
  row.names = FALSE
)

recurrent_gene_table <- de_half_sig_for_overlap %>%
  group_by(gene) %>%
  summarise(
    n_contrasts = n_distinct(set_id),
    n_celltypes = n_distinct(celltype),
    n_stages = n_distinct(stage),
    n_adaxial_up = sum(side_direction == "Adaxial up"),
    n_abaxial_up = sum(side_direction == "Abaxial up"),
    dominant_side = names(sort(table(side_direction), decreasing = TRUE))[1],
    direction_consistency = max(table(side_direction)) / n(),
    min_padj = min(p_val_adj, na.rm = TRUE),
    max_abs_log2FC = max(abs(avg_log2FC), na.rm = TRUE),
    contrasts = paste(sort(unique(set_id)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    recurrence_class = case_when(
      n_contrasts == 1 ~ "exclusive_one_contrast",
      n_celltypes >= 2 & n_stages >= 2 ~ "recurrent_celltype_and_stage",
      n_celltypes >= 2 ~ "recurrent_multiple_celltypes",
      n_stages >= 2 ~ "recurrent_multiple_stages",
      TRUE ~ "recurrent_other"
    )
  ) %>%
  arrange(desc(n_contrasts), desc(direction_consistency), min_padj)

## Add annotation if the annotated final-candidate table is available in memory.
if (exists("candidate_half") && any(c("Functional Annotation", "Best Arabidopsis BLAST hit") %in% colnames(candidate_half))) {
  recurrent_gene_anno <- candidate_half %>%
    select(any_of(c(
      "gene", "Functional Annotation", "Best Arabidopsis BLAST hit",
      "Best P. trichocarpa BLAST hit", "Genes of Interest C=cell wall",
      "Genes of Interest S=shoot apical meristem", "Genes of Interest X=xylan"
    ))) %>%
    distinct(gene, .keep_all = TRUE)
  
  recurrent_gene_table <- recurrent_gene_table %>%
    left_join(recurrent_gene_anno, by = "gene")
}

write.csv(
  recurrent_gene_table,
  "tables/DE_adaxial_abaxial_half_split_recurrent_gene_table.csv",
  row.names = FALSE
)

p_recurrence_summary <- recurrent_gene_table %>%
  count(recurrence_class, dominant_side, name = "n_genes") %>%
  ggplot(aes(x = reorder(recurrence_class, n_genes), y = n_genes, fill = dominant_side)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Adaxial up" = "#00BFC4", "Abaxial up" = "#F8766D")) +
  theme_classic(base_size = 9, base_family = overlap_base_family) +
  labs(x = NULL, y = "Number of genes", fill = "Dominant side",
       title = "Half-split side-DEG recurrence classes")

save_overlap_pdf(
  p_recurrence_summary,
  "overlap/DE_adaxial_abaxial_half_split_recurrence_summary.pdf",
  width = 7,
  height = 4
)

## 9.12.6b UpSet-style overlap matrix for half-split DEG sets ----
overlap_membership_long <- de_half_sig_for_overlap %>%
  distinct(gene, set_id, set_id_clean)

write.csv(
  overlap_membership_long,
  "tables/DE_adaxial_abaxial_half_split_overlap_membership_long.csv",
  row.names = FALSE
)

overlap_membership_wide <- overlap_membership_long %>%
  mutate(present = 1L) %>%
  tidyr::pivot_wider(
    id_cols = gene,
    names_from = set_id_clean,
    values_from = present,
    values_fill = 0L
  )

write.csv(
  overlap_membership_wide,
  "tables/DE_adaxial_abaxial_half_split_overlap_membership_matrix.csv",
  row.names = FALSE
)

set_cols <- setdiff(colnames(overlap_membership_wide), "gene")
set_sizes <- overlap_membership_wide %>%
  summarise(across(all_of(set_cols), sum)) %>%
  tidyr::pivot_longer(everything(), names_to = "set_id_clean", values_to = "n_genes") %>%
  left_join(half_split_set_key, by = "set_id_clean") %>%
  arrange(desc(n_genes))

write.csv(
  set_sizes,
  "tables/DE_adaxial_abaxial_half_split_overlap_set_sizes.csv",
  row.names = FALSE
)

## Always write a robust ggplot fallback summary. UpSetR can draw blank pages
## with some ggplot/grid/font combinations, so treat it as optional.
intersection_summary <- overlap_membership_wide %>%
  rowwise() %>%
  mutate(
    n_sets = sum(c_across(all_of(set_cols))),
    intersection = paste(set_cols[as.logical(c_across(all_of(set_cols)))], collapse = " & ")
  ) %>%
  ungroup() %>%
  count(intersection, n_sets, name = "n_genes") %>%
  filter(n_sets > 0) %>%
  arrange(desc(n_genes), desc(n_sets)) %>%
  slice_head(n = 40)

write.csv(
  intersection_summary,
  "tables/DE_adaxial_abaxial_half_split_top_intersections.csv",
  row.names = FALSE
)

p_intersections <- intersection_summary %>%
  mutate(intersection = factor(intersection, levels = rev(intersection))) %>%
  ggplot(aes(x = n_genes, y = intersection, fill = n_sets)) +
  geom_col(width = 0.75) +
  scale_fill_gradient(low = "#BFD3E6", high = "#045A8D") +
  theme_classic(base_size = 7, base_family = overlap_base_family) +
  labs(x = "Number of genes", y = NULL, fill = "Sets",
       title = "Top half-split DEG intersections",
       subtitle = "Readable fallback/companion to UpSetR")

save_overlap_pdf(
  p_intersections,
  "overlap/DE_adaxial_abaxial_half_split_top_intersections_fallback.pdf",
  width = 9,
  height = 8
)

## Use UpSetR as an optional extra if available.
if (requireNamespace("UpSetR", quietly = TRUE) && length(set_cols) >= 2) {
  upset_input <- overlap_membership_wide %>%
    select(all_of(set_cols)) %>%
    as.data.frame()
  
  pdf("overlap/DE_adaxial_abaxial_half_split_upset.pdf",
      width = 11, height = 7, family = overlap_base_family,
      useDingbats = FALSE, paper = "special")
  print(UpSetR::upset(
    upset_input,
    nsets = min(length(set_cols), 20),
    nintersects = 40,
    order.by = "freq",
    sets.bar.color = "#595959",
    main.bar.color = "#2C7FB8",
    text.scale = c(1.2, 1.1, 1, 1, 1.1, 1)
  ))
  dev.off()
}

## 9.12.6c recurrent-gene presence matrix for gene-plot selection ----
recurrent_gene_presence <- de_half_sig_for_overlap %>%
  semi_join(recurrent_gene_table %>% filter(n_contrasts >= 2), by = "gene") %>%
  mutate(
    direction_symbol = ifelse(side_direction == "Adaxial up", "A", "B"),
    signed_log2FC = avg_log2FC
  ) %>%
  select(gene, celltype, stage, side_direction, p_val_adj, avg_log2FC, set_id) %>%
  arrange(gene, celltype, stage, side_direction)

write.csv(
  recurrent_gene_presence,
  "tables/DE_adaxial_abaxial_half_split_recurrent_gene_presence_long.csv",
  row.names = FALSE
)

top_recurrent_for_matrix <- recurrent_gene_table %>%
  filter(n_contrasts >= 2) %>%
  slice_head(n = 60) %>%
  pull(gene)

presence_matrix_plot_df <- de_half_sig_for_overlap %>%
  filter(gene %in% top_recurrent_for_matrix) %>%
  mutate(
    gene = factor(gene, levels = rev(top_recurrent_for_matrix)),
    set_label = paste(celltype, stage, side_direction, sep = " | "),
    set_label = factor(set_label, levels = unique(set_label))
  )

if (nrow(presence_matrix_plot_df) > 0) {
  p_presence_matrix <- ggplot(presence_matrix_plot_df, aes(x = set_label, y = gene)) +
    geom_point(aes(size = abs(avg_log2FC), fill = side_direction),
               shape = 21, colour = "grey25", linewidth = 0.15) +
    scale_fill_manual(values = c("Adaxial up" = "#00BFC4", "Abaxial up" = "#F8766D")) +
    scale_size_continuous(range = c(1.2, 4.0)) +
    theme_classic(base_size = 6, base_family = overlap_base_family) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_blank()) +
    labs(fill = "Side", size = "|log2FC|",
         title = "Top recurrent half-split side-DEGs")
  
  save_overlap_pdf(
    p_presence_matrix,
    "overlap/DE_adaxial_abaxial_half_split_recurrent_gene_presence_matrix.pdf",
    width = 12,
    height = 9
  )
}

## 9.12.6d Mfuzz / fallback dynamic clustering of all side-biased DEGs ----
## Input is a stage pseudobulk expression matrix for all significant half-split DEGs.
## If Mfuzz is not installed, use k-means on z-scored stage profiles as a fallback.
mfuzz_priority_genes <- c(
    "PtXaTreH.09G107100.v5.1",
    "PtXaTreH.06G108400.v5.1",
    "PtXaAlbH.05G168000.v5.1",
    "PtXaTreH.05G167000.v5.1"
)

mfuzz_candidate_genes <- de_half_sig_for_overlap %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  pull(gene) %>%
  unique() %>%
  union(mfuzz_priority_genes)

mfuzz_candidate_genes <- intersect(mfuzz_candidate_genes, rownames(sobj_cross_ann))

get_avg_expression_by_group <- function(obj, genes, group_col, assay = "SCT") {
  expr <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "data")
  )
  genes <- intersect(genes, rownames(expr))
  groups <- obj@meta.data[[group_col]]
  group_levels <- unique(as.character(groups[!is.na(groups)]))
  
  avg_mat <- sapply(group_levels, function(g) {
    cells_g <- rownames(obj@meta.data)[as.character(groups) == g]
    if (length(cells_g) == 0) return(rep(NA_real_, length(genes)))
    Matrix::rowMeans(expr[genes, cells_g, drop = FALSE])
  })
  
  avg_mat <- as.matrix(avg_mat)
  rownames(avg_mat) <- genes
  avg_mat
}

if (length(mfuzz_candidate_genes) >= 10) {
  sobj_cross_ann$mfuzz_stage <- factor(as.character(sobj_cross_ann$stage),
                                       levels = stage_order_overlap, ordered = TRUE)
  mfuzz_stage_mat <- get_avg_expression_by_group(
    sobj_cross_ann,
    genes = mfuzz_candidate_genes,
    group_col = "mfuzz_stage",
    assay = "SCT"
  )
  
  ## force biological stage column order
  mfuzz_stage_mat <- mfuzz_stage_mat[, intersect(stage_order_overlap, colnames(mfuzz_stage_mat)), drop = FALSE]
  mfuzz_stage_mat <- mfuzz_stage_mat[apply(mfuzz_stage_mat, 1, function(x) sd(x, na.rm = TRUE) > 0), , drop = FALSE]
  
  write.csv(
    mfuzz_stage_mat,
    "tables/Mfuzz_half_split_all_DEG_stage_pseudobulk_matrix.csv"
  )
  write.csv(
    mfuzz_stage_mat,
    "tables/Mfuzz_half_split_recurrent_gene_stage_pseudobulk_matrix.csv"
  )
  
  mfuzz_stage_z <- t(scale(t(mfuzz_stage_mat)))
  mfuzz_stage_z[!is.finite(mfuzz_stage_z)] <- 0
  
  write.csv(
    mfuzz_stage_z,
    "tables/Mfuzz_half_split_all_DEG_stage_zscore_matrix.csv"
  )
  write.csv(
    mfuzz_stage_z,
    "tables/Mfuzz_half_split_recurrent_gene_stage_zscore_matrix.csv"
  )
  
  n_dynamic_clusters <- min(6, max(2, floor(nrow(mfuzz_stage_z) / 10)))
  
  if (requireNamespace("Mfuzz", quietly = TRUE) && requireNamespace("Biobase", quietly = TRUE) && nrow(mfuzz_stage_z) >= 10) {
    ## Mfuzz plotting functions call Biobase::exprs() as exprs() internally,
    ## so Biobase/Mfuzz need to be attached, not only available by namespace.
    suppressPackageStartupMessages({
      library(Biobase)
      library(Mfuzz)
    })
    
    eset <- Biobase::ExpressionSet(assayData = as.matrix(mfuzz_stage_mat))
    eset <- Mfuzz::filter.NA(eset, thres = 0.25)
    eset <- Mfuzz::fill.NA(eset, mode = "mean")
    eset <- Mfuzz::filter.std(eset, min.std = 0)
    eset <- Mfuzz::standardise(eset)
    m_est <- Mfuzz::mestimate(eset)
    mfuzz_res <- Mfuzz::mfuzz(eset, c = n_dynamic_clusters, m = m_est)
    
    mfuzz_max_membership <- apply(mfuzz_res$membership, 1, max, na.rm = TRUE)
    
    mfuzz_membership <- as.data.frame(mfuzz_res$membership) %>%
      rownames_to_column("gene") %>%
      mutate(
        cluster = unname(mfuzz_res$cluster[match(gene, names(mfuzz_res$cluster))]),
        dynamic_cluster_method = "Mfuzz",
        dynamic_cluster = paste0("Mfuzz_", cluster),
        max_membership = unname(mfuzz_max_membership[match(gene, names(mfuzz_max_membership))])
      )
    
    write.csv(
      mfuzz_membership,
      "tables/Mfuzz_half_split_all_DEG_membership.csv",
      row.names = FALSE
    )
    write.csv(
      mfuzz_membership,
      "tables/Mfuzz_half_split_recurrent_gene_membership.csv",
      row.names = FALSE
    )
    
    ## Reliable saved plot: rebuild the Mfuzz stage profiles with ggplot.
    ## This avoids Mfuzz::mfuzz.plot2() drawing only to the interactive device
    ## in some RStudio/macOS sessions.
    mfuzz_expr_plot_df <- as.data.frame(Biobase::exprs(eset)) %>%
      rownames_to_column("gene") %>%
      left_join(mfuzz_membership %>% select(gene, cluster), by = "gene") %>%
      tidyr::pivot_longer(
        cols = all_of(colnames(Biobase::exprs(eset))),
        names_to = "stage",
        values_to = "z"
      ) %>%
      mutate(
        stage = factor(stage, levels = colnames(Biobase::exprs(eset)), ordered = TRUE),
        cluster = factor(cluster)
      )
    
    mfuzz_cluster_means <- mfuzz_expr_plot_df %>%
      group_by(cluster, stage) %>%
      summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop")
    
    p_mfuzz_patterns <- ggplot(mfuzz_expr_plot_df, aes(stage, z, group = gene)) +
      geom_line(alpha = 0.18, linewidth = 0.25, colour = "grey45") +
      geom_line(
        data = mfuzz_cluster_means,
        aes(stage, mean_z, group = cluster),
        inherit.aes = FALSE,
        linewidth = 0.9,
        colour = "#D95F02"
      ) +
      facet_wrap(~ cluster, ncol = 3) +
      theme_classic(base_size = 8, base_family = overlap_base_family) +
      labs(
        x = "Stage",
        y = "Standardized expression",
        title = "Mfuzz stage-expression modules for all half-split DEGs",
        subtitle = "Grey lines: genes; orange line: module mean"
      )
    
    save_overlap_pdf(
      p_mfuzz_patterns,
      "mfuzz/Mfuzz_half_split_all_DEG_stage_clusters.pdf",
      width = 9,
      height = 6
    )
    save_overlap_pdf(
      p_mfuzz_patterns,
      "mfuzz/Mfuzz_half_split_recurrent_gene_stage_clusters.pdf",
      width = 9,
      height = 6
    )
    
    ## Optional native Mfuzz plot. This is kept as a secondary export because it
    ## may draw to the interactive device in some sessions.
    native_file <- "mfuzz/Mfuzz_half_split_all_DEG_stage_clusters_native.pdf"
    native_tmp <- paste0(native_file, ".tmp")
    if (file.exists(native_tmp)) unlink(native_tmp)
    native_dev_before <- grDevices::dev.cur()
    tryCatch(
      {
        grDevices::pdf(native_tmp, width = 9, height = 7, family = "sans",
                       useDingbats = FALSE, paper = "special")
        Mfuzz::mfuzz.plot2(
          eset,
          cl = mfuzz_res,
          mfrow = c(2, ceiling(n_dynamic_clusters / 2)),
          time.labels = colnames(Biobase::exprs(eset))
        )
      },
      error = function(e) {
        message("Native Mfuzz PDF export skipped: ", conditionMessage(e))
      },
      finally = {
        if (grDevices::dev.cur() != native_dev_before) grDevices::dev.off()
      }
    )
    if (file.exists(native_tmp) && file.info(native_tmp)$size > 1000) {
      if (file.exists(native_file)) unlink(native_file)
      file.rename(native_tmp, native_file)
      legacy_native_file <- "mfuzz/Mfuzz_half_split_recurrent_gene_stage_clusters_native.pdf"
      if (file.exists(legacy_native_file)) unlink(legacy_native_file)
      file.copy(native_file, legacy_native_file)
    } else if (file.exists(native_tmp)) {
      unlink(native_tmp)
    }
  } else {
    set.seed(123)
    km <- kmeans(mfuzz_stage_z, centers = n_dynamic_clusters, nstart = 50)
    
    kmeans_membership <- tibble(
      gene = rownames(mfuzz_stage_z),
      cluster = km$cluster,
      dynamic_cluster_method = "kmeans_fallback",
      dynamic_cluster = paste0("Kmeans_", cluster),
      max_membership = NA_real_
    ) %>%
      left_join(recurrent_gene_table, by = "gene")
    
    write.csv(
      kmeans_membership,
      "tables/Kmeans_fallback_half_split_all_DEG_stage_clusters.csv",
      row.names = FALSE
    )
    write.csv(
      kmeans_membership,
      "tables/Kmeans_fallback_half_split_recurrent_gene_stage_clusters.csv",
      row.names = FALSE
    )
    
    kmeans_plot_df <- as.data.frame(mfuzz_stage_z) %>%
      rownames_to_column("gene") %>%
      left_join(kmeans_membership %>% select(gene, cluster), by = "gene") %>%
      tidyr::pivot_longer(cols = all_of(colnames(mfuzz_stage_z)),
                          names_to = "stage", values_to = "z") %>%
      mutate(stage = factor(stage, levels = stage_order_overlap, ordered = TRUE))
    
    p_kmeans_patterns <- ggplot(kmeans_plot_df, aes(stage, z, group = gene)) +
      geom_line(alpha = 0.18, linewidth = 0.25, colour = "grey40") +
      stat_summary(aes(group = cluster), fun = mean, geom = "line",
                   linewidth = 0.8, colour = "#D95F02") +
      facet_wrap(~ cluster, ncol = 3) +
      theme_classic(base_size = 8, base_family = "Arial") +
      labs(x = "Stage", y = "Z-scored pseudobulk expression",
           title = "Fallback k-means stage-expression modules for all half-split DEGs",
           subtitle = "Install Mfuzz + Biobase for fuzzy c-means clustering")
    
    save_overlap_pdf(
      p_kmeans_patterns,
      "mfuzz/Kmeans_fallback_half_split_all_DEG_stage_patterns.pdf",
      width = 9,
      height = 6
    )
    save_overlap_pdf(
      p_kmeans_patterns,
      "mfuzz/Kmeans_fallback_half_split_recurrent_gene_stage_patterns.pdf",
      width = 9,
      height = 6
    )
  }
  
  ## Gene-level dynamic pattern annotation for downstream candidate tables.
  ## This adds interpretable labels such as early_peak, L4_L5_peak, late_peak,
  ## or stable_broad, regardless of whether Mfuzz or fallback k-means was used.
  dynamic_membership_for_annotation <- if (exists("mfuzz_membership")) {
    mfuzz_membership %>%
      select(gene, dynamic_cluster_method, dynamic_cluster, max_membership)
  } else if (exists("kmeans_membership")) {
    kmeans_membership %>%
      select(gene, dynamic_cluster_method, dynamic_cluster, max_membership)
  } else {
    tibble(gene = rownames(mfuzz_stage_z))
  }
  
  dynamic_gene_pattern_table <- as.data.frame(mfuzz_stage_z) %>%
    rownames_to_column("gene") %>%
    tidyr::pivot_longer(cols = all_of(colnames(mfuzz_stage_z)),
                        names_to = "stage", values_to = "z") %>%
    group_by(gene) %>%
    summarise(
      dynamic_peak_stage = stage[which.max(z)],
      dynamic_min_stage = stage[which.min(z)],
      dynamic_z_range = max(z, na.rm = TRUE) - min(z, na.rm = TRUE),
      dynamic_max_z = max(z, na.rm = TRUE),
      dynamic_min_z = min(z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      dynamic_pattern_class = case_when(
        dynamic_z_range < 0.75 ~ "stable_broad",
        dynamic_peak_stage %in% c("L1", "L2") ~ "early_peak",
        dynamic_peak_stage %in% c("L4", "L5") ~ "L4_L5_peak",
        dynamic_peak_stage %in% c("L15") ~ "late_peak",
        TRUE ~ "other_dynamic"
      )
    ) %>%
    left_join(dynamic_membership_for_annotation, by = "gene")
  
  ## Add raw mean expression across stages for easier candidate triage.
  dynamic_raw_summary <- as.data.frame(mfuzz_stage_mat) %>%
    rownames_to_column("gene") %>%
    mutate(
      dynamic_mean_expression = rowMeans(across(all_of(colnames(mfuzz_stage_mat))), na.rm = TRUE),
      dynamic_max_expression = apply(across(all_of(colnames(mfuzz_stage_mat))), 1, max, na.rm = TRUE)
    ) %>%
    select(gene, dynamic_mean_expression, dynamic_max_expression)
  
  dynamic_gene_pattern_table <- dynamic_gene_pattern_table %>%
    left_join(dynamic_raw_summary, by = "gene")
  
  write.csv(
    dynamic_gene_pattern_table,
    "tables/DE_adaxial_abaxial_half_split_dynamic_gene_pattern_table.csv",
    row.names = FALSE
  )
}

### 9.13 Venn plot: half split vs quarter split DEGs ####

library(dplyr)
library(ggVennDiagram)

deg_half <- read.csv("DE_adaxial_abaxial_significant_candidates.csv")
deg_quarter <- read.csv("DE_adaxial_abaxial_quarter_significant_candidates.csv")

# gene-level Venn
venn_gene <- list(half = unique(deg_half$gene),quater = unique(deg_quarter$gene))

p1 <- ggVennDiagram(venn_gene, label_alpha = 0,set_size = 3 ) +
  scale_fill_gradient(low = "white", high = "#55b6ff") +
  theme_void() +
  ggtitle("DEG overlap: gene level") +
  theme(plot.margin = margin(20, 20, 20, 80))

ggsave("DE_adaxial_abaxial_half_vs_quarter_gene_venn.pdf", p1, width = 8, height = 5)

# gene + celltype + stage Venn
deg_half$key <- paste(deg_half$gene, deg_half$celltype, deg_half$stage, sep = "|")
deg_quarter$key <- paste(deg_quarter$gene, deg_quarter$celltype, deg_quarter$stage, sep = "|")

venn_key <- list(
  half = unique(deg_half$key),
  quarter = unique(deg_quarter$key)
)

p2 <- ggVennDiagram(venn_key, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#55b6ff") +
  theme_void() +
  ggtitle("DEG overlap: gene-celltype-stage")

ggsave("DE_adaxial_abaxial_half_vs_quarter_gene_celltype_stage_venn.pdf", p2, width = 8, height = 5)

# export overlap tables
write.csv(data.frame(gene = intersect(venn_gene$Half_split, venn_gene$Quarter_split)),
          "tables/DE_adaxial_abaxial_half_vs_quarter_gene_overlap.csv", row.names = FALSE)

write.csv(data.frame(key = intersect(venn_key$Half_split, venn_key$Quarter_split)),
          "tables/DE_adaxial_abaxial_half_vs_quarter_gene_celltype_stage_overlap.csv", row.names = FALSE)


### 9.14 Add gene annotation to final side candidates ####

library(dplyr)
library(readr)

anno <- read.delim(file.path(ATLAS_WORK_ROOT, "genome_files", 'pop_717_hap1_2_func_anno_w_og_matrix_v6.txt'), check.names = FALSE)

final_side <- read.csv("tables/DE_adaxial_abaxial_final_spatial_side_candidates.csv")

recurrence_annotation <- if (exists("recurrent_gene_table")) {
  recurrent_gene_table
} else if (file.exists("tables/DE_adaxial_abaxial_half_split_recurrent_gene_table.csv")) {
  read.csv("tables/DE_adaxial_abaxial_half_split_recurrent_gene_table.csv")
} else {
  NULL
}

if (!is.null(recurrence_annotation) && nrow(recurrence_annotation) > 0) {
  recurrence_annotation <- recurrence_annotation %>%
    select(any_of(c(
      "gene", "n_contrasts", "n_celltypes", "n_stages",
      "n_adaxial_up", "n_abaxial_up", "dominant_side",
      "direction_consistency", "min_padj", "max_abs_log2FC",
      "recurrence_class", "contrasts"
    ))) %>%
    distinct(gene, .keep_all = TRUE) %>%
    rename(
      recurrence_n_contrasts = n_contrasts,
      recurrence_n_celltypes = n_celltypes,
      recurrence_n_stages = n_stages,
      recurrence_n_adaxial_up = n_adaxial_up,
      recurrence_n_abaxial_up = n_abaxial_up,
      recurrence_dominant_side = dominant_side,
      recurrence_direction_consistency = direction_consistency,
      recurrence_min_padj = min_padj,
      recurrence_max_abs_log2FC = max_abs_log2FC,
      recurrence_class = recurrence_class,
      recurrence_contrasts = contrasts
    )
}

dynamic_annotation <- if (exists("dynamic_gene_pattern_table")) {
  dynamic_gene_pattern_table
} else if (file.exists("tables/DE_adaxial_abaxial_half_split_dynamic_gene_pattern_table.csv")) {
  read.csv("tables/DE_adaxial_abaxial_half_split_dynamic_gene_pattern_table.csv")
} else {
  NULL
}

if (!is.null(dynamic_annotation) && nrow(dynamic_annotation) > 0) {
  dynamic_annotation <- dynamic_annotation %>%
    select(any_of(c(
      "gene", "dynamic_cluster_method", "dynamic_cluster",
      "max_membership", "dynamic_peak_stage", "dynamic_min_stage",
      "dynamic_z_range", "dynamic_max_z", "dynamic_min_z",
      "dynamic_pattern_class", "dynamic_mean_expression",
      "dynamic_max_expression"
    ))) %>%
    distinct(gene, .keep_all = TRUE)
}

anno2 <- anno %>%
  rename(gene_id = `Gene ID`) %>%
  mutate(gene = paste0(gene_id, ".v5.1")) %>%
  select(gene, gene_id, `Functional Annotation`, `Genome Location`,
         `Best Arabidopsis BLAST hit`, `Best P. trichocarpa BLAST hit`,
         `PFAM matches`, `Interpro assigned GO terms`,
         `717 recipricoal syntelog`, `P. trichocarpa syntelog`, Orthogroup_id,
         `iTAK Shiu Protein Kinase Families`,
         `iTAK Transcription Factor (TF)/Transcription Regulator (TR) Classification`,
         starts_with("Genes of Interest"))

final_side_anno <- final_side

if (!is.null(recurrence_annotation) && nrow(recurrence_annotation) > 0) {
  final_side_anno <- final_side_anno %>%
    left_join(recurrence_annotation, by = "gene")
}

if (!is.null(dynamic_annotation) && nrow(dynamic_annotation) > 0) {
  final_side_anno <- final_side_anno %>%
    left_join(dynamic_annotation, by = "gene")
}

final_side_anno <- final_side_anno %>%
  left_join(anno2, by = "gene")

write.csv(final_side_anno,
          "tables/DE_adaxial_abaxial_final_spatial_side_candidates_annotated.csv",
          row.names = FALSE)

top_plot_candidates_annotated <- final_side_anno %>%
  arrange(p_adj_side_score, p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 50)

write.csv(
  top_plot_candidates_annotated,
  "tables/DE_adaxial_abaxial_top50_plot_candidates.csv",
  row.names = FALSE
)

write.csv(
  top_plot_candidates_annotated,
  "tables/DE_adaxial_abaxial_top50_plot_candidates_annotated.csv",
  row.names = FALSE
)


qs::qsave(sobj_cross_ann, "saved_obj/sobj_petiole_cross_res0.35_final_v1_2026.7.24.qs")

  },
  mfuzz = function() {
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(tidyr);library(Matrix);library(edgeR)})
b<-file.path(ATLAS_WORK_ROOT)
src<-file.path(b,'summary/Fig5_petiole/two_tissue_development_v2');o<-file.path(src,'Mfuzz_side_comparison_v2');dir.create(o,recursive=TRUE,showWarnings=FALSE)
s<-qread(file.path(b,'saved_obj/sobj_petiole_cross_res0.35_final_v1_2026.7.24.qs'));md<-s@meta.data;md$spot<-rownames(md);s[['Spatial']]<-JoinLayers(s[['Spatial']]);mat<-GetAssayData(s,assay='Spatial',layer='counts')
cand <- read.csv(file.path(b, 'tables/DE_adaxial_abaxial_final_spatial_side_candidates_annotated.csv')) |>
 filter(method == 'split_half', celltype %in% c('Epidermis','Cortex'))
write.csv(cand, file.path(src, 'downstream_epidermis_cortex_prioritized_303.csv'), row.names=FALSE)
output<-list();qc<-list()
for(ct in c('Epidermis','Cortex')){
 mm<-md |> filter(celltypes==ct,!is.na(adaxial_abaxial));mm$key<-paste(mm$parent_library,mm$adaxial_abaxial,sep='__')
 sm<-mm |> group_by(key,parent_library,stage,adaxial_abaxial) |> summarise(n_spots=n(),.groups='drop')
 z<-sparseMatrix(i=match(mm$spot,colnames(mat)),j=match(mm$key,sm$key),x=1,dims=c(ncol(mat),nrow(sm)));pb<-mat%*%z;colnames(pb)<-sm$key
 y<-calcNormFactors(DGEList(pb));cp<-cpm(y);sm$library_size<-y$samples$lib.size;sm$TMM_factor<-y$samples$norm.factors;sm$celltype<-ct;qc[[ct]]<-sm
 genes<-unique(cand$gene[cand$celltype==ct]);stopifnot(all(genes %in% rownames(cp)))
 long<-as.data.frame(as.table(cp[genes,,drop=FALSE]));names(long)<-c('gene','key','CPM');long<-long |> mutate(gene=as.character(gene),key=as.character(key)) |> left_join(sm,by='key')
 output[[ct]]<-list(counts=pb,TMM_factors=y$samples$norm.factors,CPM=cp[genes,,drop=FALSE],meta=sm,long=long,developmental_genes=character(),side_DEGs=unique(cand$gene[cand$celltype==ct]),expressed_universe=rownames(cp)[rowSums(cp>=1)>=2])
}
saveRDS(output,file.path(o,'side_resolved_TMM_profiles.rds'));write.csv(bind_rows(qc),file.path(o,'library_side_QC.csv'),row.names=FALSE);write.csv(bind_rows(lapply(output,`[[`,'long')),file.path(o,'side_resolved_TMM_profiles.csv'),row.names=FALSE)

suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(Mfuzz);library(Biobase);library(mclust);library(patchwork)})
b<-file.path(ATLAS_WORK_ROOT);src<-file.path(b,'summary/Fig5_petiole/two_tissue_development_v2');o<-file.path(src,'Mfuzz_side_comparison_v2');dat<-readRDS(file.path(o,'side_resolved_TMM_profiles.rds'));sts<-c('L1','L2','L4','L5','L15');sides<-c('adaxial','abaxial')
export<-function(p,n,w=180,h=120){for(ex in c('pdf','svg','png')){tmp<-file.path(tempdir(),paste0(n,'.',ex));ggsave(tmp,p,width=w,height=h,units='mm',device=switch(ex,pdf=cairo_pdf,svg=svglite::svglite,png=ragg::agg_png),dpi=300,bg='white',limitsize=FALSE);stopifnot(file.info(tmp)$size>100);file.copy(tmp,file.path(o,basename(tmp)),overwrite=TRUE)}}
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold'),legend.position='top'))
sqdist<-function(x,c) pmax(outer(rowSums(x*x),rowSums(c*c),'+')-2*tcrossprod(x,c),1e-12)
assign_new<-function(x,centers,m){d<-sqdist(x,centers)^(-1/(m-1));d/rowSums(d)}
a<-read.csv(file.path(b,'genome_files/go2gene_poplar.csv'),header=FALSE);names(a)<-c('ID','gene');a$gene<-ifelse(grepl('\\.v5.1$',a$gene),a$gene,paste0(a$gene,'.v5.1'));nm<-read.csv(file.path(b,'genome_files/go2term.csv'));names(nm)<-c('ID','Description','Ontology');a<-inner_join(a,nm,by='ID') |> filter(Ontology=='biological_process') |> distinct(ID,gene,Description)
GO<-function(genes,universe){u<-intersect(unique(universe),a$gene);g<-intersect(unique(genes),u);if(length(g)<3)return(NULL);aa<-a |> filter(gene %in% u);tt<-split(aa$gene,aa$ID);tt<-tt[lengths(tt)>=10 & lengths(tt)<=500];if(!length(tt))return(NULL);cnt<-vapply(tt,function(t)sum(g %in% t),integer(1));p<-phyper(cnt-1,lengths(tt),length(u)-lengths(tt),length(g),lower.tail=FALSE);data.frame(ID=names(tt),Description=nm$Description[match(names(tt),nm$ID)],Count=cnt,group_annotated=length(g),background_annotated=length(u),term_size=lengths(tt),pvalue=p,FDR=p.adjust(p,'BH'),geneID=vapply(tt,function(t)paste(intersect(g,t),collapse='/'),character(1)))}
models<-list();metrics<-list();centroids<-list();gores<-list();members<-list();restarts<-list();sidecmp<-list()
for(ct in names(dat)){
 d<-dat[[ct]];av<-d$long |> group_by(gene,stage,adaxial_abaxial) |> summarise(CPM=mean(CPM),.groups='drop');write.csv(av,file.path(o,paste0(ct,'_side_means.csv')),row.names=FALSE)
 wide<-av |> mutate(col=paste(adaxial_abaxial,stage,sep='_')) |> select(gene,col,CPM) |> pivot_wider(names_from=col,values_from=CPM);cols<-unlist(lapply(sides,function(s)paste(s,sts,sep='_')));raw<-as.matrix(wide[,cols]);rownames(raw)<-wide$gene;lograw<-log2(raw+1)
 ad<-lograw[,paste0('adaxial_',sts)];ab<-lograw[,paste0('abaxial_',sts)];cmp<-data.frame(gene=rownames(raw),celltype=ct,shape_r=vapply(seq_len(nrow(raw)),function(i)if(sd(ad[i,])>0 && sd(ab[i,])>0)cor(ad[i,],ab[i,]) else NA_real_,numeric(1)),adaxial_range=apply(ad,1,function(x)diff(range(x))),abaxial_range=apply(ab,1,function(x)diff(range(x))),peak_adaxial=sts[max.col(ad)],peak_abaxial=sts[max.col(ab)],mean_log2_abaxial_minus_adaxial=rowMeans(ab-ad),bias_range=apply(ab-ad,1,function(x)diff(range(x))),prioritized_side_DEG=rownames(raw)%in%d$side_DEGs,developmental_gene=rownames(raw)%in%d$developmental_genes)
 cmp$both_dynamic<-cmp$adaxial_range>=1 & cmp$abaxial_range>=1;sidecmp[[ct]]<-cmp
 # Pool standardized gene-side trajectories to obtain a shared pattern dictionary for comparing the same gene across sides.
 dg<-intersect(d$side_DEGs,rownames(raw));joint<-lograw[dg,,drop=FALSE];keep<-apply(raw[dg,,drop=FALSE],1,max)>=10 & apply(joint,1,function(z)diff(range(z)))>=1;joint<-joint[keep,,drop=FALSE]
 setups<-list(Side_DEG_joint=list(x=joint,already=FALSE,meta=data.frame(id=rownames(joint),gene=rownames(joint),side='paired'),ks=2:6))
 for(mode in names(setups)){
  z<-setups[[mode]];xx<-if(z$already)z$x else t(scale(t(z$x)));stopifnot(all(is.finite(xx)));es<-ExpressionSet(xx);m<-mestimate(es);key<-paste(ct,mode,sep='_');message(key,' n=',nrow(xx),' m=',round(m,3));model<-list(x=xx,meta=z$meta,m=m,raw=if(mode=='Side_DEG_joint')raw[rownames(xx),,drop=FALSE] else NULL,fits=list())
  for(k in z$ks){
   ff<-lapply(1:15,function(j){set.seed(9162026+k*100+j);mfuzz(es,c=k,m=m,iter.max=500)});obj<-vapply(ff,function(f)sum(f$membership^m*sqdist(xx,f$centers)),numeric(1));f<-ff[[which.min(obj)]]
   ord<-if(mode=='Side_DEG_joint')order(max.col((f$centers[,1:5]+f$centers[,6:10])/2),rowMeans(f$centers[,6:10]-f$centers[,1:5])) else order(max.col(f$centers),f$centers[,1]);f$centers<-f$centers[ord,,drop=FALSE];f$membership<-f$membership[,ord,drop=FALSE];f$cluster<-max.col(f$membership);rownames(f$centers)<-paste0('G',1:k);colnames(f$membership)<-rownames(f$centers);ari<-vapply(ff,function(g)adjustedRandIndex(f$cluster,g$cluster),numeric(1));empt<-vapply(ff,function(g)sum(colSums(g$membership>=.5)==0),integer(1))
   # Sensitivity to removing 20% of genes; never interpreted as biological replication.
   boot<-sapply(1:8,function(j){set.seed(717000+k*100+j);ii<-sample(seq_len(nrow(xx)),floor(.8*nrow(xx)));g<-mfuzz(ExpressionSet(xx[ii,,drop=FALSE]),c=k,m=m,iter.max=500);pred<-max.col(assign_new(xx,g$centers,m));adjustedRandIndex(f$cluster,pred)})
   minD<-min(dist(f$centers));pc<-mean(rowSums(f$membership^2));xb<-min(obj)/(nrow(xx)*minD^2);met<-data.frame(celltype=ct,analysis=mode,k=k,n_profiles=nrow(xx),m=m,objective_per_profile=min(obj)/nrow(xx),min_centroid_distance=minD,Xie_Beni=xb,normalized_partition_coefficient=(pc-1/k)/(1-1/k),mean_restart_ARI=mean(ari),min_restart_ARI=min(ari),mean_subsample_ARI=mean(boot),min_subsample_ARI=min(boot),empty_core_runs=sum(empt>0),min_core_profiles=min(colSums(f$membership>=.5)))
   metrics[[paste(key,k)]]<-met;restarts[[paste(key,k)]]<-data.frame(celltype=ct,analysis=mode,k=k,run=1:15,objective=obj,restart_ARI=ari,empty_cores=empt)
   mem<-cbind(z$meta,data.frame(cluster=paste0('G',f$cluster),membership=apply(f$membership,1,max),celltype=ct,analysis=mode,k=k));members[[paste(key,k)]]<-mem
   cen<-as.data.frame(f$centers);cen$cluster<-rownames(f$centers);cen<-cen |> pivot_longer(-cluster,names_to='position',values_to='z') |> mutate(celltype=ct,analysis=mode,k=k,side=if(mode=='Side_DEG_joint')sub('_.*','',position) else 'combined',stage=if(mode=='Side_DEG_joint')sub('^[^_]+_','',position) else position);centroids[[paste(key,k)]]<-cen
   univ<-unique(z$meta$gene);gr<-bind_rows(lapply(paste0('G',1:k),function(cl){g<-GO(mem$gene[mem$cluster==cl & mem$membership>=.5],univ);if(!is.null(g))mutate(g,cluster=cl)}));if(nrow(gr)){gr$FDR_across_clusters<-p.adjust(gr$pvalue,'BH');gr$celltype<-ct;gr$analysis<-mode;gr$k<-k;gores[[paste(key,k)]]<-gr}
   model$fits[[as.character(k)]]<-f
   message('  k=',k,' ARI=',round(mean(ari),3),' subsample=',round(mean(boot),3),' XB=',round(xb,3),' coremin=',min(colSums(f$membership>=.5)))
  };models[[key]]<-model;saveRDS(models,file.path(o,'cluster_comparison_models.rds'))
 }
}
met<-bind_rows(metrics);cen<-bind_rows(centroids);mem<-bind_rows(members);gr<-bind_rows(gores);cmp<-bind_rows(sidecmp)
write.csv(met,file.path(o,'cluster_number_diagnostics.csv'),row.names=FALSE);write.csv(bind_rows(restarts),file.path(o,'cluster_restart_diagnostics.csv'),row.names=FALSE);write.csv(cen,file.path(o,'all_cluster_centroids.csv'),row.names=FALSE);write.csv(mem,file.path(o,'all_cluster_memberships.csv'),row.names=FALSE);write.csv(gr,file.path(o,'all_cluster_GO_tests.csv'),row.names=FALSE);write.csv(gr |> filter(FDR_across_clusters<.05,Count>=3),file.path(o,'significant_cluster_GO.csv'),row.names=FALSE);write.csv(cmp,file.path(o,'same_gene_side_comparison.csv'),row.names=FALSE)
for(mode in unique(cen$analysis)){
 cc<-cen |> filter(analysis==mode) |> mutate(stage=factor(stage,levels=sts),solution=paste0('k = ',k));p<-ggplot(cc,aes(stage,z,group=interaction(cluster,side),color=cluster,linetype=side))+geom_line(linewidth=.65)+facet_grid(celltype~solution)+labs(x='Leaf position',y='Cluster centroid (standardized expression)',title=gsub('_',' ',mode))+theme(legend.position='bottom');export(p,paste0('Cluster_number_',mode,'_v2'),240,100)
 dd<-met |> filter(analysis==mode) |> select(celltype,k,Xie_Beni,mean_restart_ARI,mean_subsample_ARI,min_centroid_distance) |> pivot_longer(-c(celltype,k),names_to='diagnostic',values_to='value');p<-ggplot(dd,aes(k,value,color=celltype))+geom_line()+geom_point()+facet_wrap(~diagnostic,scales='free_y',ncol=2)+labs(x='Number of clusters',y=NULL,title=gsub('_',' ',mode));export(p,paste0('Diagnostics_',mode,'_v2'),160,115)
}
top<-gr |> filter(FDR_across_clusters<.05,Count>=3) |> group_by(celltype,analysis,k,cluster) |> arrange(FDR_across_clusters) |> slice_head(n=4) |> ungroup();write.csv(top,file.path(o,'top_GO_by_cluster_number.csv'),row.names=FALSE)

suppressPackageStartupMessages({library(dplyr)})
b<-file.path(ATLAS_WORK_ROOT);o<-file.path(b,'summary/Fig5_petiole/two_tissue_development_v2/Mfuzz_side_comparison_v2');dat<-readRDS(file.path(o,'side_resolved_TMM_profiles.rds'));mem<-read.csv(file.path(o,'all_cluster_memberships.csv'));a<-read.csv(file.path(b,'genome_files/go2gene_poplar.csv'),header=FALSE);names(a)<-c('ID','gene');a$gene<-ifelse(grepl('\\.v5.1$',a$gene),a$gene,paste0(a$gene,'.v5.1'));nm<-read.csv(file.path(b,'genome_files/go2term.csv'));names(nm)<-c('ID','Description','Ontology');a<-inner_join(a,nm,by='ID') |> filter(Ontology=='biological_process') |> distinct(ID,gene,Description)
res<-list()
for(ct in names(dat)){
 u<-intersect(dat[[ct]]$expressed_universe,a$gene);aa<-a |> filter(gene %in% u);tt<-split(aa$gene,aa$ID);tt<-tt[lengths(tt)>=10 & lengths(tt)<=500]
 for(k in 2:6){rr<-list();for(cl in paste0('G',1:k)){
 gg<-mem |> filter(celltype==ct,analysis=='Side_DEG_joint',.data$k==!!k,cluster==cl,membership>=.5);g<-intersect(gg$gene,u);cnt<-vapply(tt,function(t)sum(g %in% t),integer(1));p<-phyper(cnt-1,lengths(tt),length(u)-lengths(tt),length(g),lower.tail=FALSE);rr[[cl]]<-data.frame(celltype=ct,k=k,cluster=cl,ID=names(tt),Description=nm$Description[match(names(tt),nm$ID)],Count=cnt,group_annotated=length(g),background_annotated=length(u),term_size=lengths(tt),pvalue=p,FDR=p.adjust(p,'BH'),geneID=vapply(tt,function(t)paste(intersect(g,t),collapse='/'),character(1)))
 };rr<-bind_rows(rr);rr$FDR_across_clusters<-p.adjust(rr$pvalue,'BH');res[[paste(ct,k)]]<-rr}
}

  }
)
atlas_dispatch(stages)
