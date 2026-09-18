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

# ======================================= #
# Bud spatial atlas pipeline (updated)   #
# ======================================= #

# ---------- local helper overrides ----------
process_spatial_integration_v5 <- function(sobj, resolution = 0.5, dims = 1:30, harmony_vars = "orig.ident") {
  sobj <- add_basic_qc(sobj)
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

rename_clusters_from_table_direct <- function(
    sobj, mapping_df,
    cluster_col = "cluster_label",
    label_col = "celltype_final",
    source_ident_col = "seurat_clusters",
    new_meta_col = "celltypes"
) {
  mapping_df <- mapping_df %>%
    mutate(
      !!cluster_col := as.character(.data[[cluster_col]]),
      !!label_col := as.character(.data[[label_col]])
    ) %>%
    filter(!is.na(.data[[cluster_col]]), !is.na(.data[[label_col]])) %>%
    distinct(.data[[cluster_col]], .keep_all = TRUE)

  old_ids <- as.character(sobj[[source_ident_col]][, 1])
  map_idx <- match(old_ids, mapping_df[[cluster_col]])
  new_ids <- mapping_df[[label_col]][map_idx]
  new_ids[is.na(new_ids)] <- old_ids[is.na(new_ids)]

  sobj[[new_meta_col]] <- new_ids
  Idents(sobj) <- sobj[[new_meta_col]][, 1]
  sobj
}

plot_known_markers_paged <- function(
    sobj, marker_genes, out_pdf, images = NULL,
    genes_per_page = 3, ncol = NULL,
    pt.size.factor = 1.8, crop = TRUE,
    width = 14, height = 10
) {
  marker_genes <- unique(marker_genes)
  marker_genes <- marker_genes[marker_genes %in% rownames(sobj)]
  if (length(marker_genes) == 0) return(invisible(NULL))

  gene_pages <- split(marker_genes, ceiling(seq_along(marker_genes) / genes_per_page))
  if (is.null(ncol)) ncol <- if (is.null(images)) genes_per_page else length(images)

  pdf(out_pdf, width = width, height = height)
  for (i in seq_along(gene_pages)) {
    gset <- gene_pages[[i]]
    p <- SpatialFeaturePlot(
      sobj, features = gset, images = images, ncol = ncol,
      crop = crop, alpha = c(0.1, 1), pt.size.factor = pt.size.factor,
      min.cutoff = "q05", max.cutoff = "q95"
    ) + plot_annotation(title = paste0("Markers (page ", i, "): ", paste(gset, collapse = ", ")))
    print(p)
  }
  dev.off()
}

run_comparecluster_from_markers <- function(marker_df, groups, universe, go_anno_bp, out_pdf, n_top = 100, title = "GO comparison") {
  gene_sets <- setNames(
    lapply(groups, function(ct) {
      marker_df %>%
        filter(cluster == ct, p_val_adj < 0.05, avg_log2FC > 0.5) %>%
        slice_max(order_by = avg_log2FC, n = n_top, with_ties = FALSE) %>%
        pull(gene) %>%
        intersect(universe)
    }),
    groups
  )
  gene_sets <- gene_sets[sapply(gene_sets, length) > 0]
  if (length(gene_sets) == 0) return(NULL)

  res <- compareCluster(
    geneCluster = gene_sets,
    fun = "enricher",
    universe = universe,
    TERM2GENE = go_anno_bp[, c("ID", "Gene")],
    TERM2NAME = go_anno_bp[, c("ID", "Description")],
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  pdf(out_pdf, width = 11, height = 8)
  print(
    dotplot(res, showCategory = 12) + ggtitle(title) &
      theme(legend.position = "right", axis.text.y = element_text(size = 10))
  )
  dev.off()
  res
}

run_full_nmf <- function(sobj, k = 12, prefix = "BUD") {
  DefaultAssay(sobj) <- "SCT"
  nmf_features <- intersect(VariableFeatures(sobj), rownames(sobj))
  mat <- GetAssayData(sobj, assay = "SCT", layer = "data")[nmf_features, ] |> as.matrix()

  set.seed(123)
  fit <- RcppML::nmf(mat, k = k)
  W <- fit$w; H <- fit$h
  rownames(W) <- rownames(mat)
  colnames(W) <- paste0("NMF_", seq_len(k))
  rownames(H) <- paste0("NMF_", seq_len(k))
  colnames(H) <- colnames(mat)

  for (f in rownames(H)) sobj[[f]] <- H[f, colnames(sobj)]

  write.csv(W, paste0("tables/", prefix, "_NMF_gene_loadings_k", k, ".csv"))
  write.csv(H, paste0("tables/", prefix, "_NMF_cell_scores_k", k, ".csv"))
  saveRDS(fit, paste0("saved_obj/", prefix, "_nmf_fit_k", k, ".rds"))

  top_genes <- purrr::map_dfr(colnames(W), function(f) {
    tibble(gene = rownames(W), loading = W[, f], factor = f) %>%
      arrange(desc(loading)) %>% slice_head(n = 50)
  })
  write.csv(top_genes, paste0("tables/", prefix, "_NMF_top50_genes.csv"), row.names = FALSE)
  list(object = sobj, W = W, H = H, fit = fit, top_genes = top_genes)
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

# ---------- 1. Load and split bud libraries ----------
raw_list <- list(
  bud_A = load_slice("data/dormant_bud_A/", "bud_A", tissue = "axillary_bud"),
  bud_B = load_slice("data/dormant_bud_B/", "bud_B", tissue = "axillary_bud"),
  bud_C = load_slice("data/dormant_bud_C/", "bud_C", tissue = "axillary_bud"),
  bud_D = load_slice("data/dormant_bud_D/", "bud_D", tissue = "axillary_bud")
)

split_params <- tibble::tribble(
  ~base_name, ~eps, ~minPts,
  "bud_A", 260, 8,
  "bud_B", 260, 8,
  "bud_C", 260, 8,
  "bud_D", 260, 8
)

pdf("QC/dbscan_section_detection_bud.pdf", width = 12, height = 6)
for (i in seq_len(nrow(split_params))) {
  nm <- split_params$base_name[i]
  print(plot_dbscan_check(raw_list[[nm]], eps = split_params$eps[i], minPts = split_params$minPts[i], title = nm))
}
dev.off()

split_res <- split_libraries_from_params(raw_list, split_params)
write.csv(split_res$centroids, "tables/bud_split_section_centroids.csv", row.names = FALSE)
section_list <- split_res$objects

# ---------- 2. Merge ----------
sobj_bud <- merge(section_list[[1]], y = section_list[-1], add.cell.ids = names(section_list), project = "Bud_split")
qsave(sobj_bud, "saved_obj/sobj_bud_split_raw.qs")

# ---------- 3. QC + integration ----------
sobj_bud <- process_spatial_integration_v5(sobj_bud, resolution = 0.5, dims = 1:30, harmony_vars = "orig.ident")
qsave(sobj_bud, "saved_obj/sobj_bud_split_harmony_res0.5.qs")

plot1 <- VlnPlot_scCustom(sobj_bud, features = c("nCount_Spatial", "nFeature_Spatial"), group.by = "orig.ident", plot_median = TRUE) + NoLegend()
plot2 <- SpatialFeaturePlot(sobj_bud, features = c("nCount_Spatial", "nFeature_Spatial"), crop = FALSE, ncol = 4)
pdf("QC/QC_bud_split_violin.pdf", width = 11, height = 7); print(plot1); dev.off()
pdf("QC/QC_bud_split_spatial.pdf", width = 14, height = 10); print(plot2); dev.off()
write.csv(describeBy(sobj_bud@meta.data, group = sobj_bud@meta.data$orig.ident, mat = TRUE), "QC/QC_stats_bud_split.csv", row.names = FALSE)

# ---------- 4. Resolution sweep ----------
resolutions <- c(0.3, 0.4, 0.5, 0.6)
images_use <- names(sobj_bud@images)[seq_len(min(4, length(names(sobj_bud@images))))]

pdf("clustering/bud_resolution_optimization.pdf", width = 12, height = 12)
for (res in resolutions) {
  sobj_tmp <- FindClusters(sobj_bud, resolution = res)
  
  for (img in images_use) {
    p <- SpatialDimPlot(
      sobj_tmp,
      cells.highlight = CellsByIdentities(sobj_tmp),
      facet.highlight = TRUE,
      images = img,
      ncol = 5,
      pt.size.factor = 4,
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

### final resolution: 0.5 #################
sobj_bud <- FindClusters(sobj_bud, resolution = 0.5)


p1 <- DimPlot_scCustom(sobj_bud, reduction = "umap", group.by = "orig.ident", alpha = 0.6,ggplot_default_colors = TRUE) +
  ggtitle("By section")
p2 <- DimPlot(sobj_bud, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("By cluster")
p3 <- SpatialDimPlot(sobj_bud, images = names(sobj_bud@images)[c(1:3)],group.by = "seurat_clusters", crop = FALSE, ncol = 3)

pdf("clustering/bud_res0.5.pdf", width = 16, height = 8)
print(p1 + p2)
print(p3)
dev.off()


# ---------- 5. Marker scoring + provisional annotation ----------
mkr_list <- read.csv("updated_markerlist_poplar_2.8.26.csv")
bud_marker_heatmap_obj <- sobj_bud
bud_final_marker_obj_file <- "saved_obj/sobj_bud_res0.5_final_v2_2026.7.24.qs"
if (file.exists(bud_final_marker_obj_file)) {
  bud_marker_heatmap_obj <- qread(bud_final_marker_obj_file)
  DefaultAssay(bud_marker_heatmap_obj) <- "SCT"
  bud_marker_group_col <- intersect(c("celltypes_v2", "celltypes", "seurat_clusters"), colnames(bud_marker_heatmap_obj@meta.data))[1]
  if (!is.na(bud_marker_group_col)) Idents(bud_marker_heatmap_obj) <- bud_marker_heatmap_obj@meta.data[[bud_marker_group_col]]
}
mkr_list_use <- mkr_list %>% filter(GeneID.v5 %in% rownames(bud_marker_heatmap_obj))
if ("out_table" %in% names(formals(make_marker_heatmap))) {
  make_marker_heatmap(bud_marker_heatmap_obj, mkr_list_use, "marker/marker_heatmap_bud_split.pdf",
                      out_table = "marker/marker_heatmap_bud_split_row_labels_renewed.csv")
} else {
  make_marker_heatmap(bud_marker_heatmap_obj, mkr_list_use, "marker/marker_heatmap_bud_split.pdf")
}

ann_res <- score_and_assign_celltypes(sobj_bud, mkr_list_use, assay = "SCT")
sobj_bud <- ann_res$object
write.csv(ann_res$cluster_scores, "tables/bud_cluster_marker_module_scores.csv", row.names = FALSE)
write.csv(ann_res$top_assign, "tables/bud_cluster_top_predicted_celltypes.csv", row.names = FALSE)
plot_annotation_qc(sobj_bud, "clustering/bud_predicted_annotations.pdf")

manual_map_template <- ann_res$top_assign %>%
  transmute(cluster_label = cluster_label, celltype_predicted = predicted_celltype, predicted_score, celltype_final = predicted_celltype) %>%
  arrange(cluster_label)
write_csv(manual_map_template, "tables/manual_cluster_annotation_template_bud.csv")

sobj_bud$celltypes_provisional <- sobj_bud$predicted_celltype
## Idents(sobj_bud) <- sobj_bud$celltypes_provisional
## qsave(sobj_bud, "saved_obj/sobj_bud_split_provisionalAnno.qs")

# ---------- 6. Apply manual annotation if available ----------
if (file.exists("tables/manual_cluster_annotation_template_bud.csv")) {
  manual_map <- read_csv("tables/manual_cluster_annotation_template_bud.csv", show_col_types = FALSE)
  sobj_bud_ann <- rename_clusters_from_table_direct(sobj_bud, manual_map, source_ident_col = "seurat_clusters", new_meta_col = "celltypes")
} else {
  sobj_bud_ann <- sobj_bud
  sobj_bud_ann$celltypes <- sobj_bud_ann$celltypes_provisional
  Idents(sobj_bud_ann) <- sobj_bud_ann$celltypes
}


## save annotated object ####### 
qsave(sobj_bud_ann, "saved_obj/sobj_bud_split_res0.5_annotated.qs") 
## update: refined celltype annotation

sobj_bud_ann <- atlas_annotations(sobj_bud_ann, "bud")
qs::qsave(sobj_bud_ann, "saved_obj/sobj_bud_res0.5_final_v2_2026.7.24.qs")

  },
  markers = function() {
sobj_bud_ann <- qs::qread(atlas_input("saved_obj/sobj_bud_res0.5_final_v2_2026.7.24.qs"))
Seurat::Idents(sobj_bud_ann) <- sobj_bud_ann$celltypes
# ---------- 8. Marker discovery after final annotation ----------
Idents(sobj_bud_ann) <- sobj_bud_ann$celltypes_v2
bud_markers <- FindAllMarkers(
  sobj_bud_ann,
  assay = "SCT",
  recorrect_umi = FALSE,
  test.use = "wilcox",
  only.pos = TRUE,
  logfc.threshold = 0.5,
  min.pct = 0.1
)
write.csv(bud_markers, "marker/BUD_all_markers_by_celltype_v2_2026_6.csv", row.names = FALSE)

bud_markers_clean <- bud_markers %>% filter(grepl("^PtXa", gene))
write.csv(bud_markers_clean, "marker/BUD_all_markers_by_celltype_PtXaOnly.csv", row.names = FALSE)

top_bud_markers <- bud_markers_clean %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

pdf("marker/BUD_top_markers_dotplot.pdf", width = 8, height = 9)
print(
  DotPlot_scCustom(
    sobj_bud_ann,
    features = unique(top_bud_markers$gene),
    flip_axes = TRUE,
    scale.by = "size",
    dot.min = 0,
    dot.scale = 6,
    x_lab_rotate = TRUE
  ) + theme(axis.text.y = element_text(size = 8))
)
dev.off()


  }
)
atlas_dispatch(stages)
