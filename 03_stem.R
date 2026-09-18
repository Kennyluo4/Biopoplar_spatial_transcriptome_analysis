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
dir.create("section_split", showWarnings = FALSE, recursive = TRUE)
dir.create("trajectory", showWarnings = FALSE, recursive = TRUE)

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
  library(ComplexHeatmap)
  library(psych)
})

# ======================================= #
# Stem spatial atlas pipeline (updated)  #
# Follows bud/SAM split-pipeline workflow #
# ======================================= #

# --------------------------------------- #
# 1. Helper functions ####
# --------------------------------------- #

merge_obj_list <- function(obj_list, project) {
  obj_list <- obj_list[!vapply(obj_list, is.null, logical(1))]
  if (length(obj_list) == 0) stop("No objects to merge.")
  if (length(obj_list) == 1) return(obj_list[[1]])
  merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list), project = project)
}

get_assay_mat <- function(sobj, assay = "SCT", layer = "data") {
  tryCatch(GetAssayData(sobj, assay = assay, layer = layer),
           error = function(e) GetAssayData(sobj, assay = assay, slot = layer))
}

process_spatial_integration_v5 <- function(sobj, resolution = 0.5, dims = 1:30, harmony_vars = "orig.ident") {
  sobj <- add_basic_qc(sobj)
  sobj <- SCTransform(sobj, assay = "Spatial", verbose = FALSE)
  sobj <- RunPCA(sobj, assay = "SCT", verbose = FALSE)
  sobj <- RunHarmony(object = sobj, group.by.vars = harmony_vars, reduction = "pca",
                     reduction.save = "harmony", plot_convergence = FALSE)
  sobj <- RunUMAP(sobj, reduction = "harmony", dims = dims)
  sobj <- FindNeighbors(sobj, reduction = "harmony", dims = dims)
  sobj <- FindClusters(sobj, resolution = resolution)
  sobj
}

rename_clusters_from_table_direct <- function(sobj, mapping_df, cluster_col = "cluster_label",
                                              label_col = "celltype_final",
                                              source_ident_col = "seurat_clusters",
                                              new_meta_col = "celltypes") {
  mapping_df <- mapping_df %>%
    mutate(!!cluster_col := as.character(.data[[cluster_col]]),
           !!label_col := as.character(.data[[label_col]])) %>%
    filter(!is.na(.data[[cluster_col]]), !is.na(.data[[label_col]])) %>%
    distinct(.data[[cluster_col]], .keep_all = TRUE)
  old_ids <- as.character(sobj[[source_ident_col]][, 1])
  new_ids <- mapping_df[[label_col]][match(old_ids, mapping_df[[cluster_col]])]
  new_ids[is.na(new_ids)] <- old_ids[is.na(new_ids)]
  sobj[[new_meta_col]] <- new_ids
  Idents(sobj) <- sobj[[new_meta_col]][, 1]
  sobj
}

plot_resolution_sweep <- function(sobj, resolutions, prefix, images_use = NULL,
                                  ncol = 4, pt.size.factor = 4, width = 14, height = 12) {
  if (is.null(images_use)) images_use <- names(sobj@images)[seq_len(min(4, length(names(sobj@images))))]
  pdf(paste0("clustering/", prefix, "_resolution_optimization.pdf"), width = width, height = height)
  for (res in resolutions) {
    sobj_tmp <- FindClusters(sobj, resolution = res)
    for (img in images_use) {
      p <- SpatialDimPlot(sobj_tmp, cells.highlight = CellsByIdentities(sobj_tmp), facet.highlight = TRUE,
                          images = img, ncol = ncol, pt.size.factor = pt.size.factor,
                          alpha = 0.8, crop = FALSE, stroke = 0) +
        plot_annotation(title = paste0(prefix, " | resolution = ", res, " | ", img))
      print(p)
    }
    print(DimPlot(sobj_tmp, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
            ggtitle(paste0(prefix, " UMAP | resolution = ", res)))
  }
  dev.off()
}

plot_final_clusters <- function(sobj, prefix, images_use = NULL, ncol = 4) {
  if (is.null(images_use)) images_use <- names(sobj@images)[seq_len(min(4, length(names(sobj@images))))]
  p1 <- DimPlot_scCustom(sobj, reduction = "umap", group.by = "orig.ident", alpha = 0.6,
                         ggplot_default_colors = TRUE) + ggtitle("By section")
  p2 <- DimPlot(sobj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("By cluster")
  p3 <- SpatialDimPlot(sobj, images = images_use, group.by = "seurat_clusters", crop = FALSE,
                       ncol = ncol, pt.size.factor = 1.8)
  pdf(paste0("clustering/", prefix, "_final_clusters.pdf"), width = 16, height = 10)
  print(p1 + p2)
  print(p3)
  dev.off()
  pdf(paste0("clustering/", prefix, "_final_clusters_split_by_cluster.pdf"), width = 14, height = 14)
  for (img in images_use) {
    print(SpatialDimPlot(sobj, cells.highlight = CellsByIdentities(sobj), facet.highlight = TRUE,
                         images = img, ncol = ncol, pt.size.factor = 2, alpha = 0.8,
                         crop = FALSE, stroke = 0) + plot_annotation(title = paste0(prefix, " | ", img)))
  }
  dev.off()
}

load_marker_list_for_scoring <- function() {
  if (file.exists("updated_markerlist_poplar_2.8.26.csv")) {
    read.csv("updated_markerlist_poplar_2.8.26.csv")
  } else {
    mkr <- read.csv("marker_list_poplar_9.27.24.csv")
    mkr$GeneID.v5 <- paste0(mkr$geneID, ".v5.1")
    if (!"celltype" %in% colnames(mkr) && "cell_type" %in% colnames(mkr)) mkr$celltype <- mkr$cell_type
    mkr
  }
}

plot_known_markers_paged <- function(sobj, marker_genes, out_pdf, images = NULL,
                                     genes_per_page = 8, ncol = 4, pt.size.factor = 2,
                                     crop = FALSE, width = 16, height = 12) {
  marker_genes <- unique(marker_genes)
  marker_genes <- marker_genes[marker_genes %in% rownames(sobj)]
  if (length(marker_genes) == 0) return(invisible(NULL))
  gene_pages <- split(marker_genes, ceiling(seq_along(marker_genes) / genes_per_page))
  pdf(out_pdf, width = width, height = height)
  for (i in seq_along(gene_pages)) {
    print(SpatialFeaturePlot(sobj, features = gene_pages[[i]], images = images, ncol = ncol,
                             crop = crop, alpha = c(0.1, 1), pt.size.factor = pt.size.factor,
                             min.cutoff = "q05", max.cutoff = "q95") +
            plot_annotation(title = paste0("Known markers page ", i)))
  }
  dev.off()
}

make_manual_template_from_current_clusters <- function(sobj, ann_res, file) {
  manual_map_template <- ann_res$top_assign %>%
    transmute(cluster_label = cluster_label,
              celltype_predicted = predicted_celltype,
              predicted_score,
              celltype_final = predicted_celltype) %>%
    arrange(as.numeric(cluster_label))
  write_csv(manual_map_template, file)
  manual_map_template
}

run_comparecluster_from_markers <- function(marker_df, groups, universe, go_anno_bp,
                                            out_pdf, n_top = 100, title = "GO comparison") {
  gene_sets <- setNames(lapply(groups, function(ct) {
    marker_df %>% filter(cluster == ct, p_val_adj < 0.05, avg_log2FC > 0.5) %>%
      slice_max(order_by = avg_log2FC, n = n_top, with_ties = FALSE) %>%
      pull(gene) %>% intersect(universe)
  }), groups)
  gene_sets <- gene_sets[sapply(gene_sets, length) > 0]
  if (length(gene_sets) == 0) return(NULL)
  res <- compareCluster(geneCluster = gene_sets, fun = "enricher", universe = universe,
                        TERM2GENE = go_anno_bp[, c("ID", "Gene")],
                        TERM2NAME = go_anno_bp[, c("ID", "Description")],
                        pvalueCutoff = 0.05, pAdjustMethod = "BH")
  pdf(out_pdf, width = 12, height = 8)
  print(dotplot(res, showCategory = 12) + ggtitle(title) &
          theme(legend.position = "right", axis.text.y = element_text(size = 10)))
  dev.off()
  res
}

run_full_nmf <- function(sobj, k = 12, prefix = "STEM") {
  DefaultAssay(sobj) <- "SCT"
  nmf_features <- intersect(VariableFeatures(sobj), rownames(sobj))
  mat <- get_assay_mat(sobj, assay = "SCT", layer = "data")[nmf_features, ] |> as.matrix()
  set.seed(123)
  fit <- RcppML::nmf(mat, k = k)
  W <- fit$w; H <- fit$h
  rownames(W) <- rownames(mat); colnames(W) <- paste0("NMF_", seq_len(k))
  rownames(H) <- paste0("NMF_", seq_len(k)); colnames(H) <- colnames(mat)
  for (f in rownames(H)) sobj[[f]] <- H[f, colnames(sobj)]
  write.csv(W, paste0("tables/", prefix, "_NMF_gene_loadings_k", k, ".csv"))
  write.csv(H, paste0("tables/", prefix, "_NMF_cell_scores_k", k, ".csv"))
  saveRDS(fit, paste0("saved_obj/", prefix, "_nmf_fit_k", k, ".rds"))
  top_genes <- purrr::map_dfr(colnames(W), function(f) {
    tibble(gene = rownames(W), loading = W[, f], factor = f) %>% arrange(desc(loading)) %>% slice_head(n = 50)
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

# --------------------------------------- #
# 2. Load raw stem libraries ####
# --------------------------------------- #

raw_list <- list(
  stem_A = load_slice("data/stem2_GC_A/", "stem_A", tissue = "stem_cross"),
  stem_B = load_slice("data/stem2_GC_B/", "stem_B", tissue = "stem_cross")
)

raw_list$stem_A$section_type <- "cross"
raw_list$stem_B$section_type <- "cross"

# --------------------------------------- #
# 3. Tune section splitting per library ####
# --------------------------------------- #

# A/B are cross sections. C/D are longitudinal sections and are kept separate.
# Tune eps/minPts after inspecting section_split/dbscan_section_detection_stem_all.pdf.
split_params <- tibble::tribble(
  ~base_name, ~eps, ~minPts, ~section_type,
  "stem_A", 260, 8, "cross",
  "stem_B", 260, 8, "cross"
)

pdf("section_split/dbscan_section_detection_stem_all.pdf", width = 12, height = 8)
for (i in seq_len(nrow(split_params))) {
  nm <- split_params$base_name[i]
  print(plot_dbscan_check(raw_list[[nm]], eps = split_params$eps[i], minPts = split_params$minPts[i],
                          title = paste0(nm, " | ", split_params$section_type[i])))
}
dev.off()

# --------------------------------------- #
# 4. Split libraries into individual sections ####
# --------------------------------------- #

split_res <- split_libraries_from_params(raw_list, split_params %>% select(base_name, eps, minPts))
write.csv(split_res$centroids, "tables/stem_split_section_centroids.csv", row.names = FALSE)
section_list <- split_res$objects

for (nm in names(section_list)) {
  section_list[[nm]]$section_type <- ifelse(grepl("^stem_[AB]", nm), "cross", "longitudinal")
  section_list[[nm]]$stem_library <- sub("^(stem_[A-D]).*", "\\1", nm)
}

cross_section_list <- section_list[grepl("^stem_[AB]", names(section_list))]
long_section_list <- section_list[grepl("^stem_[CD]", names(section_list))]
# qsave(section_list, "saved_obj/stem_all_split_section_list.qs")
# qsave(cross_section_list, "saved_obj/stem_cross_AB_split_section_list.qs")
# qsave(long_section_list, "saved_obj/stem_long_CD_split_section_list.qs")

# --------------------------------------- #
# 5. Initial QC across all split sections ####
# --------------------------------------- #

sobj_stem_all_split <- merge_obj_list(section_list, project = "Stem_all_split")
sobj_stem_all_split <- add_basic_qc(sobj_stem_all_split)
# qsave(sobj_stem_all_split, "saved_obj/sobj_stem_all_split_raw.qs")

pdf("QC/QC_stem_all_split_violin.pdf", width = 12, height = 7)
print(VlnPlot_scCustom(sobj_stem_all_split, features = c("nCount_Spatial", "nFeature_Spatial"),
                       group.by = "orig.ident", plot_median = TRUE) + NoLegend())
dev.off()

pdf("QC/QC_stem_all_split_spatial.pdf", width = 16, height = 10)
print(SpatialFeaturePlot(sobj_stem_all_split, features = c("nCount_Spatial", "nFeature_Spatial"),
                         crop = FALSE, ncol = 5)) & theme(legend.position = 'right')
dev.off()
write.csv(describeBy(sobj_stem_all_split@meta.data, group = sobj_stem_all_split@meta.data$orig.ident, mat = TRUE),
          "QC/QC_stats_stem_all_split.csv", row.names = FALSE)

# --------------------------------------- #
# 6. Merge and integrate stem cross sections A/B ####
# --------------------------------------- #

sobj_stem_cross <- merge_obj_list(cross_section_list, project = "Stem_cross_AB_split")
qsave(sobj_stem_cross, "saved_obj/sobj_stem_cross_AB_split_raw.qs")

# initial resolution here only creates graph/UMAP; final resolution is chosen below after sweep
sobj_stem_cross <- process_spatial_integration_v5(sobj_stem_cross, resolution = 0.5, dims = 1:30, harmony_vars = "orig.ident")
# qsave(sobj_stem_cross, "saved_obj/sobj_stem_cross_AB_harmony_initial.qs")

# --------------------------------------- #
# 7. QC for stem cross sections A/B ####
# --------------------------------------- #

plot1 <- VlnPlot_scCustom(sobj_stem_cross, features = c("nCount_Spatial", "nFeature_Spatial"),
                          group.by = "orig.ident", plot_median = TRUE) + NoLegend()
plot2 <- SpatialFeaturePlot(sobj_stem_cross, features = c("nCount_Spatial", "nFeature_Spatial"),
                            crop = FALSE, ncol = 4)
pdf("QC/QC_stem_cross_AB_violin.pdf", width = 11, height = 7); print(plot1); dev.off()
pdf("QC/QC_stem_cross_AB_spatial.pdf", width = 14, height = 10); print(plot2); dev.off()
write.csv(describeBy(sobj_stem_cross@meta.data, group = sobj_stem_cross@meta.data$orig.ident, mat = TRUE),
          "QC/QC_stats_stem_cross_AB.csv", row.names = FALSE)

# --------------------------------------- #
# 8. Resolution sweep for stem cross sections A/B ####
# --------------------------------------- #

# Inspect this PDF first. Then update final_res_cross below.
resolutions_cross <- c(0.4, 0.5, 0.6, 0.8, 1.0)
images_cross <- names(sobj_stem_cross@images)[seq_len(min(4, length(names(sobj_stem_cross@images))))]
plot_resolution_sweep(sobj_stem_cross, resolutions_cross, prefix = "STEM_cross_AB", images_use = images_cross,
                      ncol = 3, pt.size.factor = 1.4)

### final resolution: edit after visual review #################
final_res_cross <- 0.6
sobj_stem_cross <- FindClusters(sobj_stem_cross, resolution = final_res_cross)
sobj_stem_cross$stem_cross_res0.6_cluster <- as.character(sobj_stem_cross$seurat_clusters)

## Subcluster cluster 3 to separate cambium and phloem ------------
## Reviewed temporary sweep: resolution 1.1 gives 4 subclusters.
Idents(sobj_stem_cross) <- sobj_stem_cross$stem_cross_res0.6_cluster
stem_cluster3 <- subset(sobj_stem_cross, idents = "3")
DefaultAssay(stem_cluster3) <- "SCT"
stem_cluster3_npcs <- min(20, ncol(stem_cluster3) - 1)
stem_cluster3 <- RunPCA(stem_cluster3, assay = "SCT", npcs = stem_cluster3_npcs, verbose = FALSE)
stem_cluster3 <- RunUMAP(stem_cluster3, reduction = "pca", dims = 1:stem_cluster3_npcs)
stem_cluster3 <- FindNeighbors(stem_cluster3, reduction = "pca", dims = 1:stem_cluster3_npcs)
stem_cluster3 <- FindClusters(stem_cluster3, resolution = 1.1)
stem_cluster3_sub <- as.character(stem_cluster3$seurat_clusters)
names(stem_cluster3_sub) <- colnames(stem_cluster3)
sobj_stem_cross$stem_cross_cluster3_sub_res1.1 <- NA_character_
sobj_stem_cross@meta.data[names(stem_cluster3_sub), "stem_cross_cluster3_sub_res1.1"] <- stem_cluster3_sub
stem_cross_split_clusters <- sobj_stem_cross$stem_cross_res0.6_cluster
stem_cross_split_clusters[names(stem_cluster3_sub)] <- paste0("3_", stem_cluster3_sub)
stem_cross_split_levels <- unique(stem_cross_split_clusters)
stem_cross_split_levels <- stem_cross_split_levels[order(
  as.numeric(sub("_.*", "", stem_cross_split_levels)),
  ifelse(grepl("_", stem_cross_split_levels), as.numeric(sub(".*_", "", stem_cross_split_levels)), -1)
)]
sobj_stem_cross$seurat_clusters <- factor(stem_cross_split_clusters, levels = stem_cross_split_levels)
Idents(sobj_stem_cross) <- sobj_stem_cross$seurat_clusters
qsave(sobj_stem_cross, "saved_obj/sobj_stem_cross_AB_res0.6_subcluster3_res1.1_raw.qs")
plot_final_clusters(sobj_stem_cross, prefix = paste0("STEM_cross_AB_res", final_res_cross),
                    images_use = images_cross, ncol = 4)

# --------------------------------------- #
# 9. Marker scoring + provisional annotation for stem cross sections A/B ####
# --------------------------------------- #

mkr_list <- load_marker_list_for_scoring()
mkr_list_use <- mkr_list %>% filter(GeneID.v5 %in% rownames(sobj_stem_cross))

# marker heatmap and automated marker score prediction
make_marker_heatmap(sobj_stem_cross, mkr_list_use, "marker/marker_heatmap_stem_cross_AB_split.pdf")

# check selected known markers spatially
stem_known_markers <- mkr_list_use %>%
  filter(celltype %in% c("Epidermis", "Cortex", "Pith", "Cambium", "Phloem", "Xylem", "Vascular")) %>%
  distinct(GeneID.v5) %>%
  pull(GeneID.v5) %>%
  intersect(rownames(sobj_stem_cross))

plot_known_markers_paged(
  sobj = sobj_stem_cross,
  marker_genes = stem_known_markers,
  out_pdf = "marker/STEM_cross_AB_known_markers_spatial.pdf",
  images = names(sobj_stem_cross@images),
  genes_per_page = 4,
  ncol = 4,
  pt.size.factor = 2,
  crop = FALSE,
  width = 16,
  height = 10
)


# dotplot of known markers across clusters
marker_label_df <- mkr_list_use %>%
  filter(GeneID.v5 %in% stem_known_markers) %>%
  group_by(GeneID.v5) %>%
  summarise(
    celltype = paste(unique(celltype), collapse = "/"),
    Name = paste(unique(Name), collapse = "/"),
    .groups = "drop"
  ) %>%
  mutate(plot_label = paste0(celltype, " | ", Name, " | ", GeneID.v5))

# make sure each label is unique
marker_label_df$plot_label <- make.unique(marker_label_df$plot_label)

stem_known_markers_use <- marker_label_df$GeneID.v5
label_map <- setNames(marker_label_df$plot_label, marker_label_df$GeneID.v5)

pdf("marker/STEM_cross_AB_known_markers_dotplot.pdf", width = 10, height = 14)
p <- DotPlot_scCustom(
  sobj_stem_cross,
  features = stem_known_markers_use,
  flip_axes = TRUE,
  dot.scale = 6,
  x_lab_rotate = TRUE
) +
  scale_x_discrete(labels = label_map) +
  theme(axis.text.y = element_text(size = 8))

print(p)
dev.off()

# create manual annotation template directly from current clusters
manual_map_template_cross <- tibble(
  cluster_label = sort(unique(as.character(sobj_stem_cross$seurat_clusters))),
  celltype_predicted = cluster_label,
  predicted_score = NA_real_,
  celltype_final = case_when(
    cluster_label == "0" ~ "Cortex",
    cluster_label == "1" ~ "Epidermis",
    cluster_label == "2" ~ "Pith",
    cluster_label %in% c("3_0", "3_3") ~ "Cambium",
    cluster_label %in% c("3_1", "3_2") ~ "Phloem",
    cluster_label == "4" ~ "Phloem",
    cluster_label == "5" ~ "Xylem",
    TRUE ~ cluster_label
  )
)

write_csv(
  manual_map_template_cross,
  "tables/manual_cluster_annotation_template_stem_cross_AB.csv"
)


# --------------------------------------- #
# 10. Apply manual annotation for stem cross sections A/B ####
# --------------------------------------- #

# Edit tables/manual_cluster_annotation_template_stem_cross_AB.csv manually after checking:
# 1) clustering/STEM_cross_AB_res*_final_clusters*.pdf
# 2) marker/STEM_cross_AB_known_markers_spatial_paged.pdf
# 3) marker/marker_heatmap_stem_cross_AB_split.pdf
# Suggested coarse labels for cross section: Epidermis, Cortex, Phloem/Cambium, Xylem, Pith.
# Re-run from here after editing the CSV.
manual_map_cross <- read_csv("tables/manual_cluster_annotation_template_stem_cross_AB.csv",
                             show_col_types = FALSE)

sobj_stem_cross_ann <- rename_clusters_from_table_direct(
  sobj = sobj_stem_cross,
  mapping_df = manual_map_cross,
  cluster_col = "cluster_label",
  label_col = "celltype_final",
  source_ident_col = "seurat_clusters",
  new_meta_col = "celltypes"
)

## Save the annotated object #######
qsave(sobj_stem_cross_ann, "saved_obj/sobj_stem_cross_AB_annotated_v2_2026.7.qs")


# --------------------------------------- #
# 11. Final overview of annotated stem cross sections ####
# --------------------------------------- #

pdf("clustering/STEM_cross_AB_final_annotation_dimplot.pdf", width = 9, height = 8)
print(DimPlot(sobj_stem_cross_ann, reduction = "umap", group.by = "celltypes", label = TRUE, repel = TRUE) +
        ggtitle("Stem cross-section annotated UMAP"))
dev.off()

pdf("clustering/STEM_cross_AB_final_annotation_spatialDimPlot.pdf", width = 18, height = 9)
print(SpatialDimPlot(sobj_stem_cross_ann, group.by = "celltypes", crop = FALSE, ncol = 4, pt.size.factor = 1.6) +
        ggtitle("Stem cross-section annotated spatial map"))
dev.off()

plot_cluster_proportions(sobj_stem_cross_ann, order_vec = names(sobj_stem_cross_ann@images),
                         cluster_col = "celltypes",
                         out_pdf = "clustering/STEM_cross_AB_celltype_proportions_by_section.pdf")


sobj_stem_cross_ann <- atlas_annotations(sobj_stem_cross_ann, "stem")
qs::qsave(sobj_stem_cross_ann, "saved_obj/sobj_stem_cross_AB_res0.6_final_v2_Seurat_RCTD_2026.7.26.qs")

  },
  markers = function() {
sobj_stem_cross_ann <- qs::qread(atlas_input("saved_obj/sobj_stem_cross_AB_res0.6_final_v2_Seurat_RCTD_2026.7.26.qs"))
Seurat::Idents(sobj_stem_cross_ann) <- sobj_stem_cross_ann$celltypes
# 12. Marker discovery after final stem cross-section annotation ####
# --------------------------------------- #

Idents(sobj_stem_cross_ann) <- sobj_stem_cross_ann$celltypes
stem_markers <- FindAllMarkers(sobj_stem_cross_ann, assay = "SCT", recorrect_umi = FALSE,
                               test.use = "wilcox", only.pos = TRUE,
                               logfc.threshold = 0.5, min.pct = 0.1)
write.csv(stem_markers, "marker/STEM_cross_AB_all_markers_by_celltype.csv", row.names = FALSE)

stem_markers_clean <- stem_markers %>% filter(grepl("^PtXa", gene))
write.csv(stem_markers_clean, "marker/STEM_cross_AB_all_markers_by_celltype_PtXaOnly.csv", row.names = FALSE)

top_stem_markers <- stem_markers_clean %>% group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>% ungroup()

pdf("marker/STEM_cross_AB_top_markers_dotplot.pdf", width = 8, height = 9)
print(DotPlot_scCustom(sobj_stem_cross_ann, features = unique(top_stem_markers$gene),
                       flip_axes = TRUE, scale.by = "size", dot.min = 0,
                       dot.scale = 6, x_lab_rotate = TRUE) +
        theme(axis.text.y = element_text(size = 8)))
dev.off()

pdf("marker/STEM_cross_AB_top_markers_spatial.pdf", width = 16, height = 12)
print(SpatialFeaturePlot(sobj_stem_cross_ann, features = unique(top_stem_markers$gene)[1:min(16, length(unique(top_stem_markers$gene)))],
                         images = names(sobj_stem_cross_ann@images)[1], pt.size.factor = 4,
                         crop = T, ncol = 4, min.cutoff = "q05", max.cutoff = "q95")) & theme(legend.position = 'right')
dev.off()


  }
)
atlas_dispatch(stages)
