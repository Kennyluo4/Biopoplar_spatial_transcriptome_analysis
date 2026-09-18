# Submission edition: run from this directory or set POPLAR_CODE_ROOT.
.code_root <- Sys.getenv("POPLAR_CODE_ROOT", unset = "")
if (!nzchar(.code_root)) {
  .script <- grep("^--file=", commandArgs(), value = TRUE)
  .code_root <- if (length(.script)) dirname(normalizePath(sub("^--file=", "", .script[1]))) else getwd()
}
source(file.path(.code_root, "R", "submission_setup.R"))

stages <- list(
  preprocess = function() {

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(harmony)
  library(scCustomize)
  library(qs)
  library(ComplexHeatmap)
  library(psych)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(circlize)
  library(magick)
  library(png)
  library(jsonlite)
})
walk(c("QC", "clustering", "marker", "saved_obj", "tables"), dir.create, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------- #
# 0. Helper functions ####################
# --------------------------------------- #

load_slice <- function(data_dir, slice_name, tissue_type = NULL) {
  obj <- Load10X_Spatial(data.dir = data_dir, slice = slice_name)
  obj$orig.ident <- slice_name
  obj$slice_id <- slice_name
  obj$tissue.type <- ifelse(is.null(tissue_type), slice_name, tissue_type)
  obj
}

process_spatial_integration <- function(sobj, resolution = 0.25, dims = 1:30) {
  sobj <- PercentageFeatureSet(sobj, pattern = "^ChrMt", col.name = "percent_mt")
  sobj <- PercentageFeatureSet(sobj, pattern = "^ChrPt", col.name = "percent_chlp")

  sobj <- SCTransform(sobj, assay = "Spatial", verbose = FALSE)
  sobj <- RunPCA(sobj, assay = "SCT", verbose = FALSE)
  sobj <- RunHarmony(
    object = sobj,
    group.by.vars = "orig.ident",
    assay.use = "SCT",
    reduction.save = "harmony",
    plot_convergence = FALSE
  )
  sobj <- RunUMAP(sobj, reduction = "harmony", dims = dims)
  sobj <- FindNeighbors(sobj, reduction = "harmony", dims = dims)
  sobj <- FindClusters(sobj, resolution = resolution)
  sobj
}

plot_cluster_proportions <- function(sobj, order_vec, cluster_col = NULL, out_pdf) {
  plot_data <- if (is.null(cluster_col)) {
    sobj@meta.data %>% mutate(cluster = as.character(Idents(sobj)))
  } else {
    sobj@meta.data %>% mutate(cluster = as.character(.data[[cluster_col]]))
  }

  plot_data <- plot_data %>%
    group_by(orig.ident, cluster) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(orig.ident) %>%
    mutate(proportion = n_cells / sum(n_cells))

  plot_data$orig.ident <- factor(plot_data$orig.ident, levels = order_vec)

  p <- ggplot(plot_data, aes(x = orig.ident, y = proportion, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = cluster), position = position_stack(vjust = 0.5), size = 3, color = "white") +
    theme_classic() +
    labs(x = NULL, y = "Proportion", fill = ifelse(is.null(cluster_col), "Ident", cluster_col)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(out_pdf, p, width = 11, height = 6)
}

make_marker_heatmap <- function(sobj, marker_df, out_pdf, group_col = "celltypes_v2") {
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

  md <- sobj@meta.data
  md <- md[, !grepl("^score_|^tmpScore_", colnames(md)), drop = FALSE]
  sobj@meta.data <- md

  marker_sets <- marker_df %>%
    distinct(celltype, GeneID.v5) %>%
    filter(GeneID.v5 %in% rownames(sobj)) %>%
    group_by(celltype) %>%
    summarise(genes = list(unique(GeneID.v5)), .groups = "drop")

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

  list(object = sobj, cluster_scores = cluster_scores, top_assign = top_assign, score_long = score_long)
}

plot_annotation_qc <- function(sobj, out_pdf) {
  p1 <- DimPlot(sobj, reduction = "umap", group.by = "orig.ident", alpha = 0.7) + ggtitle("By section")
  p2 <- DimPlot(sobj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("By cluster")
  p3 <- SpatialDimPlot(sobj, group.by = "seurat_clusters", crop = FALSE, ncol = 3)
  p4 <- SpatialDimPlot(sobj, group.by = "predicted_celltype", crop = FALSE, ncol = 3)
  pdf(out_pdf, width = 16, height = 10)
  print(p1 + p2)
  print(p3)
  print(p4)
  dev.off()
}

plot_known_markers <- function(
    sobj,
    marker_genes,
    out_pdf,
    images = NULL,
    genes_per_page = 3,
    ncol = NULL,
    pt.size.factor = 2,
    crop = TRUE,
    width = 14,
    height = 10
) {
  marker_genes <- unique(marker_genes)
  marker_genes <- marker_genes[marker_genes %in% rownames(sobj)]
  if (length(marker_genes) == 0) return(invisible(NULL))
  
  gene_pages <- split(marker_genes, ceiling(seq_along(marker_genes) / genes_per_page))
  
  if (is.null(ncol)) {
    ncol <- if (is.null(images)) genes_per_page else length(images)
  }
  
  pdf(out_pdf, width = width, height = height)
  
  for (i in seq_along(gene_pages)) {
    gset <- gene_pages[[i]]
    
    p <- SpatialFeaturePlot(
      sobj,
      features = gset,
      images = images,
      ncol = ncol,
      crop = crop,
      alpha = c(0.1, 1),
      pt.size.factor = pt.size.factor,
      min.cutoff = "q05",
      max.cutoff = "q95"
    ) +
      plot_annotation(
        title = paste0("Markers (page ", i, "): ", paste(gset, collapse = ", "))
      )
    
    print(p)
  }
  
  dev.off()
}

rename_clusters_from_table <- function(sobj, mapping_df, cluster_col = "cluster_label", label_col = "celltype_final") {
  mapping_df <- mapping_df %>%
    dplyr::mutate(
      !!cluster_col := as.character(.data[[cluster_col]]),
      !!label_col := as.character(.data[[label_col]])
    ) %>%
    dplyr::filter(!is.na(.data[[cluster_col]]), !is.na(.data[[label_col]])) %>%
    dplyr::distinct(.data[[cluster_col]], .keep_all = TRUE)
  
  cl_to_type <- stats::setNames(as.list(mapping_df[[label_col]]), mapping_df[[cluster_col]])
  
  sobj <- RenameIdents(sobj, !!!cl_to_type)
  sobj$celltypes <- Idents(sobj)
  sobj
}

resize_and_inject <- function(seurat_obj, image_map, ratio = 0.7) {
  for (slice in names(image_map)) {
    if (!slice %in% names(seurat_obj@images)) next
    base <- image_map[[slice]]
    orig_img <- file.path(base, "spatial/tissue_hires_image.png")
    new_img  <- file.path(base, "spatial/tissue_hires_image_lighter.png")
    img <- image_read(orig_img)
    image_write(image_scale(img, paste0(round(image_info(img)$width * ratio), "x")), path = new_img, format = "png")
    seurat_obj@images[[slice]]@image <- readPNG(new_img)
    seurat_obj@images[[slice]]@scale.factors$lowres <- fromJSON(file.path(base, "spatial/scalefactors_json.json"))$tissue_hires_scalef * ratio
    message("Updated image: ", slice)
  }
  seurat_obj
}

z_score <- function(data_matrix) {
  mean_values <- rowMeans(data_matrix)
  sd_values <- apply(data_matrix, 1, sd)
  (data_matrix - mean_values) / sd_values
}

cor_zscore <- function(df1, df2, variable_gene_pct = 0.9) {
  var1 <- apply(df1, 1, mad)
  var2 <- apply(df2, 1, mad)
  cutoff1 <- quantile(var1, probs = variable_gene_pct)
  cutoff2 <- quantile(var2, probs = variable_gene_pct)
  df1.filter <- df1[var1 > cutoff1, , drop = FALSE]
  df2.filter <- df2[var2 > cutoff2, , drop = FALSE]
  common_genes <- intersect(rownames(df1.filter), rownames(df2.filter))
  z1 <- z_score(df1[common_genes, , drop = FALSE])
  z2 <- z_score(df2[common_genes, , drop = FALSE])
  t(cor(as.matrix(z1), as.matrix(z2), method = "spearman"))
}

# --------------------------------------- #
# 1. Load longitudinal libraries #########
# --------------------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(harmony)
  library(scCustomize)
  library(qs)
  library(psych)
})
dir.create("QC", showWarnings = FALSE, recursive = TRUE)
dir.create("clustering", showWarnings = FALSE, recursive = TRUE)
dir.create("marker", showWarnings = FALSE, recursive = TRUE)
dir.create("saved_obj", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

load_slice_with_meta <- function(data_dir, slice_name, stage, image_quality, batch) {
  obj <- Load10X_Spatial(data.dir = data_dir, slice = slice_name)
  obj$orig.ident <- slice_name
  obj$slice_id <- slice_name
  obj$stage <- stage
  obj$image_quality <- image_quality
  obj$batch <- batch
  obj$tissue_type <- "petiole_longitudinal"
  obj
}

long_specs <- tribble(
  ~data_dir,                 ~slice_name,                  ~stage, ~image_quality,  ~batch,
  "data/petiole1_A_sub2/",   "petiole_l1_longitudinal1",  "L1",   "good",          "1st",
  "data/petiole1_B_sub1/",   "petiole_l1_longitudinal2",  "L1",   "good",          "1st",
  "data/petiole1_C_sub1/",   "petiole_l1_longitudinal3",  "L1",   "good",          "1st",
  "data/petiole2_B_sub1/",   "petiole_l2_longitudinal",   "L2",   "good",          "1st",
  "data/petiole1_D_sub1/",   "petiole_l4_longitudinal1",  "L4",   "good",          "1st",
  "data/petiole2_A_sub1/",   "petiole_l4_longitudinal2",  "L4",   "good",          "1st",
  "data/petiole2_C/",        "petiole_l5_longitudinal1",  "L5",   "damaged_image", "1st",
  "data/petiole3_B_L5long/", "petiole_l5_longitudinal2",  "L5",   "good",          "2nd",
  "data/petiole2_D/",        "petiole_l15_longitudinal",  "L15",  "low_quality",   "1st"
)

long_list <- pmap(
  long_specs,
  ~ load_slice_with_meta(..1, ..2, ..3, ..4, ..5)
)
names(long_list) <- long_specs$slice_name

write.csv(long_specs, "tables/petiole_longitudinal_library_specs.csv", row.names = FALSE)


# --------------------------------------- #
# 2. Merge objects #######################
# --------------------------------------- #

sobj_long <- merge(
  long_list[[1]],
  y = long_list[-1],
  add.cell.ids = names(long_list),
  project = "Petiole_Longitudinal"
)

## update metadata
sobj_long$stage <- factor(sobj_long$stage, levels = c("L1", "L2", "L4", "L5", "L15"))
sobj_long$image_quality <- factor(sobj_long$image_quality,
                                  levels = c("good", "damaged_image", "low_quality"))
sobj_long$batch <- factor(sobj_long$batch, levels = c("1st", "2nd"))

# qsave(sobj_long, "saved_obj/sobj_longitudinal_raw_petiole_crossStyle.qs")

rm(long_list)
gc()

# --------------------------------------- #
# 3. QC + integration ####################
# --------------------------------------- #

process_spatial_integration <- function(sobj, resolution = 0.25, dims = 1:30,
                                        harmony_vars = "orig.ident") {
  sobj <- PercentageFeatureSet(sobj, pattern = "^ChrMt", col.name = "percent_mt")
  sobj <- PercentageFeatureSet(sobj, pattern = "^ChrPt", col.name = "percent_chlp")
  
  sobj <- SCTransform(sobj, assay = "Spatial", verbose = FALSE)
  sobj <- RunPCA(sobj, assay = "SCT", verbose = FALSE)
  
  sobj <- RunHarmony(
    object = sobj,
    group.by.vars = harmony_vars,
    reduction = "pca",
    reduction.save = "harmony",
    plot_convergence = FALSE
  )
  
  sobj <- RunUMAP(sobj, reduction = "harmony", dims = dims)
  sobj <- FindNeighbors(sobj, reduction = "harmony", dims = dims)
  sobj <- FindClusters(sobj, resolution = resolution)
  
  sobj
}

## use orig.ident for primary integration
## batch is kept in metadata for QC/diagnostics because only one L5 library is in batch 2
sobj_long <- process_spatial_integration(
  sobj_long,
  resolution = 0.25,
  dims = 1:30,
  harmony_vars = "orig.ident"
)

qsave(sobj_long, "saved_obj/sobj_longitudinal_harmony_res0.15_petiole_crossStyle.qs")
# --------------------------------------- #
# 4. QC plots and batch diagnostics ######
# --------------------------------------- #

plot1 <- VlnPlot_scCustom(
  sobj_long,
  features = c("nCount_Spatial", "nFeature_Spatial"),
  group.by = "orig.ident",
  plot_median = TRUE
) + NoLegend()

plot2 <- SpatialFeaturePlot(
  sobj_long,
  features = c("nCount_Spatial", "nFeature_Spatial"),
  crop = FALSE,
  ncol = 5,
  pt.size.factor = 1.5
) 
# & theme(legend.position = "right")

pdf("QC/QC_violin_petiole_longitudinal.pdf", width = 11, height = 7)
print(plot1)
dev.off()

pdf("QC/QC_spatial_petiole_longitudinal.pdf", width = 16, height = 14)
print(plot2)
dev.off()

sum_stats <- describeBy(
  sobj_long@meta.data,
  group = sobj_long@meta.data$orig.ident,
  mat = TRUE
)
write.csv(sum_stats, "QC/QC_stats_beforeFilter_petiole_longitudinal.csv", row.names = FALSE)

## UMAP diagnostics: section / stage / batch / image quality
p_umap1 <- DimPlot(sobj_long, reduction = "umap", group.by = "orig.ident", label = FALSE) +
  ggtitle("UMAP by section")
p_umap2 <- DimPlot(sobj_long, reduction = "umap", group.by = "stage", label = TRUE) +
  ggtitle("UMAP by stage")
p_umap3 <- DimPlot(sobj_long, reduction = "umap", group.by = "batch", label = TRUE) +
  ggtitle("UMAP by batch")
p_umap4 <- DimPlot(sobj_long, reduction = "umap", group.by = "image_quality", label = TRUE) +
  ggtitle("UMAP by image quality")

pdf("QC/UMAP_batch_stage_quality_petiole_longitudinal.pdf", width = 14, height = 10)
print((p_umap1 | p_umap2) / (p_umap3 | p_umap4))
dev.off()

## optional: inspect only the suspicious L5 section from batch 2
pdf("QC/L5_batch2_check_petiole_longitudinal.pdf", width = 10, height = 6)
print(
  SpatialDimPlot(
    sobj_long,
    group.by = "seurat_clusters",
    images = "petiole_l5_longitudinal2",
    pt.size.factor = 1.5,
    crop = FALSE
  ) +
    plot_annotation(title = "petiole_l5_longitudinal2 (batch 2)")
)
dev.off()

# --------------------------------------- #
# 5. Resolution sweep ####################
# --------------------------------------- #

resolutions <- c(0.15, 0.2, 2.5, 0.3, 0.4, 0.5)

## choose representative sections:
## - good L1
## - good L4
## - good L5 from batch 2
## - low-quality L15
images_use <- c(
  "petiole_l1_longitudinal2",
  "petiole_l4_longitudinal1",
  "petiole_l5_longitudinal2",
  "petiole_l15_longitudinal"
)

pdf("clustering/longitudinal_resolution_optimization_splitClusters_petiole.pdf", width = 14, height = 12)

for (res in resolutions) {
  sobj_tmp <- FindClusters(sobj_long, resolution = res)
  
  for (img in images_use) {
    p <- SpatialDimPlot(
      sobj_tmp,
      cells.highlight = CellsByIdentities(sobj_tmp),
      facet.highlight = TRUE,
      images = img,
      ncol = 4,
      pt.size.factor = 2,
      alpha = 0.8
    ) +
      plot_annotation(
        title = paste0("Resolution = ", res, " | ", img),
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
    
    print(p)
  }
}
dev.off()

## after visual review, set final resolution
### final resolution: 0.15 #################
sobj_long <- FindClusters(sobj_long, resolution = 0.15)

p1 <- DimPlot_scCustom(sobj_long, reduction = "umap", group.by = "orig.ident", alpha = 0.6,ggplot_default_colors = TRUE) +
  ggtitle("By section")
p2 <- DimPlot(sobj_long, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("By cluster")
p3 <- SpatialDimPlot(sobj_long, group.by = "seurat_clusters", crop = FALSE, ncol = 3)

pdf("clustering/petiole_longitudinal_integration_res0.15.pdf", width = 16, height = 8)
print(p1 + p2)
print(p3)
dev.off()

qsave(sobj_long, "saved_obj/sobj_longitudinal_harmony_res0.15_petiole.qs")

# --------------------------------------- #
# 6. Marker-based annotation #############
# --------------------------------------- #

## Load the final resolution object:
# sobj_long <- qread("saved_obj/sobj_long_harmony_res0.15.qs")

mkr_list <- read.csv("updated_markerlist_poplar_2.8.26.csv")
mkr_list_use <- mkr_list %>% filter(GeneID.v5 %in% rownames(sobj_long))

make_marker_heatmap(
  sobj = sobj_long,
  marker_df = mkr_list_use,
  out_pdf = "marker/marker_heatmap_longitudinal_petiole_cross.pdf"
)

ann_res <- score_and_assign_celltypes(sobj_long, mkr_list_use, assay = "SCT")
sobj_long <- ann_res$object

write.csv(ann_res$cluster_scores, "tables/petiole_longitudinal_cluster_marker_module_scores_petiole_cross.csv", row.names = FALSE)
write.csv(ann_res$top_assign, "tables/petiole_longitudinal_cluster_top_predicted_celltypes_petiole_cross.csv", row.names = FALSE)
write.csv(ann_res$score_long, "tables/petiole_longitudinal_cluster_marker_module_scores_long_petiole_cross.csv", row.names = FALSE)

plot_annotation_qc(
  sobj_long,
  out_pdf = "clustering/petiole_longitudinal_predicted_annotations_petiole_cross.pdf"
)

marker_groups <- list(
  Cortex    = c("Cortex", "Ground meristem", "Endodermis"),
  Epidermis = c("Epidermis", "Guard cell", "Trichomes"),
  Phloem    = c("Phloem", "Phloem Mother Cell", "Companion Cell",
                "Sieve Element", "Sieve element - companion cell complex"),
  Xylem     = c("Xylem", "Vessel elements", "Ray cells"),
  Vascular  = c("Cambium", "Vascular"),
  Other     = c("Pith", "Cork", "Meristematic", "Shoot meristematic", "Proliferating")
)

images_use <- c(
  "petiole_l1_longitudinal2",
  "petiole_l4_longitudinal1",
  "petiole_l5_longitudinal2",
  "petiole_l15_longitudinal"
)

pick_markers <- function(df, celltypes, obj, n_each = 2) {
  df %>%
    filter(celltype %in% celltypes) %>%
    distinct(celltype, GeneID.v5, .keep_all = TRUE) %>%
    filter(GeneID.v5 %in% rownames(obj)) %>%
    group_by(celltype) %>%
    slice_head(n = n_each) %>%
    ungroup() %>%
    pull(GeneID.v5) %>%
    unique()
}

for (grp in names(marker_groups)) {
  genes <- pick_markers(mkr_list_use, marker_groups[[grp]], sobj_long, n_each = 2)
  if (length(genes) == 0) next
  
  gene_pages <- split(genes, ceiling(seq_along(genes) / 3))
  
  for (i in seq_along(gene_pages)) {
    p <- SpatialFeaturePlot(
      sobj_long,
      features = gene_pages[[i]],
      images = images_use,
      ncol = length(images_use),
      crop = TRUE,
      alpha = c(0.1, 1),
      pt.size.factor = 1.8,
      min.cutoff = "q05",
      max.cutoff = "q95"
    ) + plot_annotation(title = paste0(grp, " markers"))
    
    ggsave(
      paste0("marker/", grp, "_page", i, "_petiole_longitudinal.png"),
      p, width = 14, height = 8, dpi = 300
    )
  }
}


manual_map_template <- ann_res$top_assign %>%
  transmute(cluster_label, celltype_predicted = predicted_celltype, predicted_score, celltype_final = predicted_celltype) %>%
  arrange(as.numeric(cluster_label))
write_csv(manual_map_template, "tables/manual_cluster_annotation_template_petiole_longitudinal.csv")

sobj_long$celltypes_provisional <- sobj_long$predicted_celltype
Idents(sobj_long) <- sobj_long$celltypes_provisional

pdf("clustering/petiole_longitudinal_provisional_celltypes.pdf", width = 10, height = 8)
print(DimPlot(sobj_long, reduction = "umap", group.by = "celltypes_provisional", label = TRUE))
dev.off()
pdf("clustering/petiole_longitudinal_provisional_celltypes_petiole_cross.pdf", width = 16, height = 8)
print(SpatialDimPlot(sobj_long, group.by = "celltypes_provisional", crop = FALSE, ncol = 3))
dev.off()

qsave(sobj_long, "saved_obj/sobj_longitudinal_petiole_cross_provisionalAnno.qs")

# --------------------------------------- #
# 7. Manual renaming #####################
# --------------------------------------- #
# Edit tables/manual_cluster_annotation_template_petiole_longitudinal.csv after visual review.
# Then uncomment below:


# manual_map <- read_csv("tables/manual_cluster_annotation_template_petiole_longitudinal.csv", show_col_types = FALSE)
# sobj_long <- rename_clusters_from_table(sobj_long, manual_map)

# ## or rename manually
Idents(sobj_long) <- sobj_long$seurat_clusters
sobj_long_ann <- RenameIdents(sobj_long, '0' = 'Cortex1', '1'='Cortex2', '2'='Vasculature1',
                               '3'='Epidermis','4'='Xylem1', '5'='Phloem','6'='Pith', '7'='Xylem2', 
                               '8'='Vasculature2','9'='Cortex1', '10'='Cortex1')
sobj_long_ann$celltypes <- Idents(sobj_long_ann)


qsave(sobj_long_ann, "saved_obj/sobj_longitudinal_petiole_cross_annotated.qs")

sobj_long_ann

pdf("clustering/petiole_longitudinal_res0.15_annotated_Dimplot.pdf", width = 10, height = 8)
print(DimPlot(sobj_long, reduction = "umap", group.by = "celltypes_provisional", label = TRUE))
dev.off()
pdf("clustering/petiole_longitudinal_res0.15_annotated_spatialDimplot.pdf", width = 16, height = 8)
print(SpatialDimPlot(sobj_long, group.by = "celltypes_provisional", crop = FALSE, ncol = 3))
dev.off()

# --------------------------------------- #

sobj_long_ann <- atlas_annotations(sobj_long_ann, "petiole_long")
qs::qsave(sobj_long_ann, "saved_obj/sobj_petiole_longitudinal_res0.15_final_v1_2026.7.24.qs")

  },
  markers = function() {
sobj_long_ann <- qs::qread(atlas_input("saved_obj/sobj_petiole_longitudinal_res0.15_final_v1_2026.7.24.qs"))
Seurat::Idents(sobj_long_ann) <- sobj_long_ann$celltypes
# 8. Marker discovery + SVG ##############
# --------------------------------------- #

# Idents(sobj_long_ann) <- sobj_long_ann$celltypes
all_de_markers_petiole_long <- FindAllMarkers(
  sobj_long_ann,
  recorrect_umi = FALSE,
  test.use = "wilcox",
  logfc.threshold = 1,
  only.pos = TRUE,
  min.diff.pct = 0.2,
  min.pct = 0.1
)
# all_de_markers_petiole_long <- read.csv('marker/sp_all_denovo_markers_petiole_longitudinal_petiole_cross_res0.15_logfc1_20pctDif.csv')      ## old version
write.csv(all_de_markers_petiole_long, "marker/PETIOLE_LONG_all_markers_by_celltype_PtXaOnly.csv", row.names = FALSE)


  }
)
atlas_dispatch(stages)
