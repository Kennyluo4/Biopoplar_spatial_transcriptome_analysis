# Submission edition: run from this directory or set POPLAR_CODE_ROOT.
.code_root <- Sys.getenv("POPLAR_CODE_ROOT", unset = "")
if (!nzchar(.code_root)) {
  .script <- grep("^--file=", commandArgs(), value = TRUE)
  .code_root <- if (length(.script)) dirname(normalizePath(sub("^--file=", "", .script[1]))) else getwd()
}
source(file.path(.code_root, "R", "submission_setup.R"))

stages <- list(
  figures123 = function() {
# =============================================================== #
# Poplar spatial atlas figure generation script (updated)
# =============================================================== #
# Purpose:
#   Generate figure-ready plots after tissue-specific pipelines have been split.
#   This script assumes each tissue pipeline has already created annotated Seurat
#   objects, marker tables, NMF outputs, GO outputs, and trajectory objects/results.
#
# Notes:
#   1. This script does NOT redo tissue-level analysis.
#   2. Some panels require manual stitching in Illustrator/Inkscape/PowerPoint.
#   3. Update object paths in Section 0.2 to match the latest output names from
#      each tissue-specific pipeline.
#   4. For panels using already generated PDFs, comments are left as reminders.


# --------------------------------------- #
# 0. Setup ####
# --------------------------------------- #

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(qs)
  library(scCustomize)
  library(ComplexHeatmap)
  library(circlize)
  library(Matrix)     # Fig. 3: per-domain rowMeans on sparse matrices
  library(scales)     # Fig. 3: rescale() for the shared z-score colour ramp
  library(FNN)        # Fig. 3B/3D: kNN spatial smoothing for display
  library(future)     # Fig. 3.7: PrepSCTFindMarkers needs sequential plan
})
# setwd(file.path(ATLAS_WORK_ROOT))   ## for windows 

dir.create("summary", showWarnings = FALSE, recursive = TRUE)
dir.create("summary/Fig1_overview", showWarnings = FALSE, recursive = TRUE)
dir.create("summary/Fig2_annotation", showWarnings = FALSE, recursive = TRUE)
dir.create("summary/Fig3_marker_atlas", showWarnings = FALSE, recursive = TRUE)
dir.create("summary/Fig4_bud_trajectory", showWarnings = FALSE, recursive = TRUE)
dir.create("summary/Fig5_stem_trajectory", showWarnings = FALSE, recursive = TRUE)
dir.create("summary/Fig6_programs_GO", showWarnings = FALSE, recursive = TRUE)
dir.create("summary/Supplementary", showWarnings = FALSE, recursive = TRUE)
dir.create("summary/Supplementary/FigS_petiole_axis", showWarnings = FALSE, recursive = TRUE)

# ---------- 0.1 small helper functions ----------

save_pdf <- function(plot, file, width = 8, height = 6) {
  pdf(file, width = width, height = height)
  print(plot)
  dev.off()
}

safe_qread <- function(file) {
  if (!file.exists(file)) {
    message("Missing file: ", file)
    return(NULL)
  }
  qread(file)
}

safe_read_csv <- function(file) {
  if (!file.exists(file)) {
    message("Missing file: ", file)
    return(NULL)
  }
  read.csv(file, check.names = FALSE)
}

get_meta_col <- function(obj, candidates) {
  candidates <- candidates[candidates %in% colnames(obj@meta.data)]
  if (length(candidates) == 0) return(NULL)
  candidates[1]
}

remove_images_for_merge <- function(obj) {
  for (img in Images(obj)) obj[[img]] <- NULL
  obj
}

read_feature_list <- function(file) {
  if (!file.exists(file)) return(character(0))
  x <- readLines(file, warn = FALSE)
  x <- trimws(x)
  unique(x[nzchar(x) & !grepl("^#", x)])
}

read_spaceranger_feature_aliases <- function(feature_root = "data") {
  feature_files <- list.files(
    feature_root,
    pattern = "features.tsv(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(feature_files) == 0) {
    message("No Space Ranger features.tsv files found under ", feature_root, ".")
    return(tibble(feature_id = character(), feature_name = character(), feature_name_unique = character()))
  }

  purrr::map_dfr(feature_files, function(file) {
    feature_df <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(feature_df) < 2) return(NULL)
    tibble(
      feature_id = feature_df[[1]],
      feature_name = feature_df[[2]],
      feature_name_unique = make.unique(feature_df[[2]])
    )
  }) %>%
    distinct(feature_id, feature_name, feature_name_unique)
}

match_features_for_qc <- function(obj, assay, list_file, fallback_pattern, label,
                                  feature_aliases = NULL) {
  assay_features <- rownames(obj[[assay]])
  listed_features <- read_feature_list(list_file)

  if (length(listed_features) > 0) {
    listed_aliases <- character(0)
    if (!is.null(feature_aliases) && nrow(feature_aliases) > 0) {
      listed_aliases <- feature_aliases %>%
        filter(feature_id %in% listed_features | feature_name %in% listed_features | feature_name_unique %in% listed_features) %>%
        select(feature_id, feature_name, feature_name_unique) %>%
        unlist(use.names = FALSE) %>%
        unique()
    }

    matched_features <- intersect(unique(c(listed_features, listed_aliases)), assay_features)
    message(
      "Organelle QC ", label, ": read ", length(listed_features),
      " IDs from ", list_file, "; matched ", length(matched_features),
      " assay features."
    )
    if (length(matched_features) > 0) return(matched_features)
    warning("No ", label, " IDs from ", list_file, " matched assay rownames. Falling back to pattern ", fallback_pattern, ".")
  } else {
    message("Organelle QC ", label, ": ", list_file, " not found; falling back to pattern ", fallback_pattern, ".")
  }

  grep(fallback_pattern, assay_features, value = TRUE)
}

expand_feature_aliases <- function(features, feature_aliases) {
  if (length(features) == 0) return(character(0))
  aliases <- character(0)
  if (!is.null(feature_aliases) && nrow(feature_aliases) > 0) {
    aliases <- feature_aliases %>%
      filter(feature_id %in% features | feature_name %in% features | feature_name_unique %in% features) %>%
      select(feature_id, feature_name, feature_name_unique) %>%
      unlist(use.names = FALSE) %>%
      unique()
  }
  unique(c(features, aliases))
}

add_organelle_qc <- function(obj, assay = "Spatial",
                             chrpt_file = "ChrPt_gene_ids.txt",
                             chrmt_file = "ChrMt_gene_ids.txt",
                             exclude_chrmt_features = c("gene-rrnL", "gene-rrnS", "rrnL", "rrnS")) {
  if (is.null(obj)) return(NULL)
  DefaultAssay(obj) <- assay

  feature_aliases <- read_spaceranger_feature_aliases()
  chrpt_features <- match_features_for_qc(obj, assay, chrpt_file, "^ChrPt", "ChrPt", feature_aliases = feature_aliases)
  chrmt_features <- match_features_for_qc(obj, assay, chrmt_file, "^ChrMt", "ChrMt", feature_aliases = feature_aliases)
  chrmt_features_all <- chrmt_features
  excluded_chrmt_features <- intersect(expand_feature_aliases(exclude_chrmt_features, feature_aliases), chrmt_features_all)
  chrmt_features <- setdiff(chrmt_features_all, excluded_chrmt_features)
  message("Organelle QC matched features in ", assay, ": ChrPt=", length(chrpt_features), "; ChrMt=", length(chrmt_features_all))
  message("Organelle QC ChrMt metric excludes ", length(excluded_chrmt_features), " rRNA features: ", paste(excluded_chrmt_features, collapse = ", "))

  if (length(chrpt_features) > 0) {
    obj <- PercentageFeatureSet(obj, features = chrpt_features, assay = assay, col.name = "percent_chrpt")
  } else {
    warning("No features matching '^ChrPt' were found in assay ", assay, ". percent_chrpt set to NA.")
    obj$percent_chrpt <- NA_real_
  }

  if (length(chrmt_features_all) > 0) {
    obj <- PercentageFeatureSet(obj, features = chrmt_features_all, assay = assay, col.name = "percent_chrmt_all")
  } else {
    warning("No features matching '^ChrMt' were found in assay ", assay, ". percent_chrmt_all set to NA.")
    obj$percent_chrmt_all <- NA_real_
  }

  if (length(chrmt_features) > 0) {
    obj <- PercentageFeatureSet(obj, features = chrmt_features, assay = assay, col.name = "percent_chrmt")
  } else {
    warning("No ChrMt features remained after rrnL/rrnS exclusion in assay ", assay, ". percent_chrmt set to NA.")
    obj$percent_chrmt <- NA_real_
  }

  obj$percent_chrmt_rrnLS <- obj$percent_chrmt_all - obj$percent_chrmt
  obj$percent_organelle_all <- rowSums(
    cbind(obj$percent_chrpt, obj$percent_chrmt_all),
    na.rm = TRUE
  )
  obj$percent_organelle_all[is.na(obj$percent_chrpt) & is.na(obj$percent_chrmt_all)] <- NA_real_
  obj$percent_organelle <- rowSums(
    cbind(obj$percent_chrpt, obj$percent_chrmt),
    na.rm = TRUE
  )
  obj$percent_organelle[is.na(obj$percent_chrpt) & is.na(obj$percent_chrmt)] <- NA_real_
  obj
}

plot_spatial_annotation <- function(obj, tissue_name, file, group_col = NULL, ncol = 4,
                                    width = 12, height = 10, pt.size.factor = 1.6,
                                    title_size = 10, text_size = 10, legend_size = 8) {
  if (is.null(obj)) return(NULL)
  if (is.null(group_col)) group_col <- get_meta_col(obj, c("celltypes_v2", "celltypes", "celltype", "assigned_cell_type", "seurat_clusters"))
  
  p <- SpatialDimPlot(
    obj,
    group.by = group_col,
    crop = FALSE,
    ncol = ncol,
    pt.size.factor = pt.size.factor
  ) +
    plot_annotation(title = tissue_name) &
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5), legend.text = element_text(size = legend_size),
      legend.title = element_text(size = legend_size),strip.text = element_text(size = text_size)
    )
   ggsave(file, p, width = width, height = height, device = cairo_pdf)
}
plot_umap_annotation <- function(obj, tissue_name, file, group_col = NULL,
                                 width = 8, height = 7) {
  if (is.null(obj)) return(NULL)
  if (is.null(group_col)) group_col <- get_meta_col(obj, c("celltypes", "celltype", "assigned_cell_type", "seurat_clusters"))
  p <- DimPlot(obj, reduction = "umap", group.by = group_col, label = TRUE, repel = TRUE) + ggtitle(tissue_name)
  save_pdf(p, file, width = width, height = height)
}

plot_marker_dotplot <- function(obj, marker_df, file, group_col = NULL, n = 5,
                                width = 9, height = 10) {
  if (is.null(obj) || is.null(marker_df)) return(NULL)
  if (is.null(group_col)) group_col <- get_meta_col(obj, c("celltypes", "celltype", "assigned_cell_type", "seurat_clusters"))
  if (!"gene" %in% colnames(marker_df)) return(NULL)
  if (!"cluster" %in% colnames(marker_df)) marker_df$cluster <- marker_df[[group_col]]
  fc_col <- if ("avg_log2FC" %in% colnames(marker_df)) "avg_log2FC" else if ("avg_logFC" %in% colnames(marker_df)) "avg_logFC" else NULL
  if (is.null(fc_col)) return(NULL)

  top_markers <- marker_df %>%
    filter(gene %in% rownames(obj)) %>%
    group_by(cluster) %>%
    slice_max(order_by = .data[[fc_col]], n = n, with_ties = FALSE) %>%
    ungroup()

  if (nrow(top_markers) == 0) return(NULL)
  Idents(obj) <- obj[[group_col]][, 1]
  p <- DotPlot_scCustom(obj, features = unique(top_markers$gene), flip_axes = TRUE,
                        scale.by = "size", dot.min = 0, dot.scale = 6, x_lab_rotate = TRUE) +
    theme(axis.text.y = element_text(size = 8))
  save_pdf(p, file, width = width, height = height)
}

plot_spatial_features_safe <- function(obj, features, file, ncol = 3, width = 14, height = 10,
                                       crop = FALSE, pt.size.factor = 1.6) {
  if (is.null(obj)) return(NULL)
  features <- intersect(unique(features), rownames(obj))
  if (length(features) == 0) {
    message("No valid features for ", file)
    return(NULL)
  }
  p <- SpatialFeaturePlot(obj, features = features, crop = crop, ncol = ncol,
                          pt.size.factor = pt.size.factor,
                          min.cutoff = "q05", max.cutoff = "q95")
  save_pdf(p, file, width = width, height = height)
}

plot_nmf_spatial <- function(obj, prefix, factors = paste0("NMF_", 1:8), image_use = NULL,
                             file = NULL, width = 15, height = 12) {
  if (is.null(obj)) return(NULL)
  factors <- intersect(factors, colnames(obj@meta.data))
  if (length(factors) == 0) {
    message("No NMF factors found for ", prefix)
    return(NULL)
  }
  if (is.null(file)) file <- paste0("summary/Fig6_programs_GO/", prefix, "_NMF_spatial.pdf")
  p <- SpatialFeaturePlot(obj, features = factors, images = image_use, crop = FALSE,
                          ncol = 4, pt.size.factor = 1.6) & theme(legend.position = "right")
  save_pdf(p, file, width = width, height = height)
}

plot_gene_pseudotime_line <- function(cds_use, genes, file, n_bins = 20, width = 7, height = 5) {
  genes <- intersect(genes, rownames(cds_use))
  if (length(genes) == 0) return(NULL)
  pt <- pseudotime(cds_use)
  expr_mat <- as.matrix(log1p(counts(cds_use)[genes, names(pt), drop = FALSE]))

  plot_df <- as.data.frame(t(expr_mat)) %>%
    rownames_to_column("cell") %>%
    mutate(pseudotime = pt[cell], pseudotime_bin = cut(pseudotime, breaks = n_bins)) %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") %>%
    group_by(gene, pseudotime_bin) %>%
    summarize(mean_expr = mean(expr, na.rm = TRUE),
              mean_pseudotime = mean(pseudotime, na.rm = TRUE), .groups = "drop")

  p <- ggplot(plot_df, aes(mean_pseudotime, mean_expr, color = gene)) +
    geom_line(linewidth = 1) + theme_classic() +
    labs(x = "Pseudotime", y = "Mean log1p counts")
  ggsave(file, p, width = width, height = height)
}

# --------------------------------------- #
# 0.2 Load updated tissue-specific objects ####
# --------------------------------------- #
# These final objects are generated by the individual tissue pipelines.
sobj_bud_ann <- safe_qread("saved_obj/sobj_bud_res0.5_final_v2_2026.7.24.qs")
sobj_stem_cross_ann <- safe_qread("saved_obj/sobj_stem_cross_AB_res0.6_final_v2_Seurat_RCTD_2026.7.26.qs")
sobj_sam_ann <- safe_qread("saved_obj/sobj_sam_split_cleaned_res0.4_annotated_v2_2026.7.qs")
sobj_petiole_cross_ann <- safe_qread("saved_obj/sobj_petiole_cross_res0.35_final_v1_2026.7.24.qs")
sobj_petiole_long_ann <- safe_qread("saved_obj/sobj_petiole_longitudinal_res0.15_final_v1_2026.7.24.qs")

# # Fallback to older file names from the previous summary script if final objects are not found.
# if (is.null(sobj_sam_ann)) sobj_sam_ann <- safe_qread("saved_obj/sobj_sp_SAM_GC_2024_merged_harmony_res0.5_anno_2.7.26.qs")
# if (is.null(sobj_stem_cross_ann)) sobj_stem_cross_ann <- safe_qread("saved_obj/sobj_sp_poplar_stem2_merged2_harmony_integration_res0.8_101124.qs")
# if (is.null(sobj_petiole_cross_ann)) sobj_petiole_cross_ann <- safe_qread("saved_obj/sobj_cross_harmony_res0.7.qs")
# if (is.null(sobj_petiole_long_ann)) sobj_petiole_long_ann <- safe_qread("saved_obj/sobj_long_harmony_res0.5.qs")

# assign tissue labels for atlas-level plots
if (!is.null(sobj_bud_ann)) sobj_bud_ann$tissue <- "Axillary bud"
if (!is.null(sobj_stem_cross_ann)) sobj_stem_cross_ann$tissue <- "Stem"
if (!is.null(sobj_sam_ann)) sobj_sam_ann$tissue <- "SAM"
if (!is.null(sobj_petiole_cross_ann)) sobj_petiole_cross_ann$tissue <- "Petiole cross"
if (!is.null(sobj_petiole_long_ann)) sobj_petiole_long_ann$tissue <- "Petiole longitudinal"

obj_list <- list(
  bud = sobj_bud_ann,
  stem = sobj_stem_cross_ann,
  sam = sobj_sam_ann,
  petiole_cross = sobj_petiole_cross_ann,
  petiole_long = sobj_petiole_long_ann
)
obj_list <- obj_list[!vapply(obj_list, is.null, logical(1))]

# --------------------------------------- #
# 1. Atlas overview figure ####
# --------------------------------------- #

## Fig. 1A: Experimental design cartoon --------
# MANUAL PANEL: create outside R.
# Suggested content:
#   Poplar plant -> SAM, axillary bud, stem, petiole cross, petiole longitudinal
#   Fresh tissue -> cryosection -> Visium -> H&E -> sequencing -> Seurat -> atlas

## Fig. 1B: representative tissue spatial maps, not annotation ----
plot_spatial_tissue_overview <- function(obj, tissue_name, file, tissue_color, ncol = 4, width = 16, height = 9, pt.size.factor = 1.6) {
  if (is.null(obj)) return(NULL)
  obj$.fig1_tissue <- tissue_name
  
  p <- SpatialDimPlot(
    obj,
    group.by = ".fig1_tissue",
    crop = FALSE,
    ncol = ncol,
    pt.size.factor = pt.size.factor,
    cols = tissue_color
  ) +
    NoLegend() +
    plot_annotation(title = tissue_name)
  
  save_pdf(p, file, width = width, height = height)
}

fig1b_info <- tibble::tribble(
  ~obj_name, ~tissue_name, ~file_name,
  "sobj_bud_ann", "Axillary bud", "Fig1B_bud_representative_spatial.pdf",
  "sobj_sam_ann", "SAM", "Fig1B_sam_representative_spatial.pdf",
  "sobj_stem_cross_ann", "Stem", "Fig1B_stem_representative_spatial.pdf",
  "sobj_petiole_cross_ann", "Petiole cross", "Fig1B_petiole_cross_representative_spatial.pdf",
  "sobj_petiole_long_ann", "Petiole longitudinal", "Fig1B_petiole_long_representative_spatial.pdf"
)

## define a color for each tissue, keep it consistent throughout the figure
tissue_cols <- c(
  "Axillary bud" = "#F8766D",
  "Petiole cross" = "#00B0F4",
  "Petiole longitudinal" = "#00BF7D",
  "SAM" = "#E76BF3",
  "Stem" = "#A3A500"
)
tissue_cols_merged_petiole <- c(
  "Axillary bud" = unname(tissue_cols["Axillary bud"]),
  "Petiole" = "#00B0F4",
  "SAM" = unname(tissue_cols["SAM"]),
  "Stem" = unname(tissue_cols["Stem"])
)

purrr::pwalk(fig1b_info, function(obj_name, tissue_name, file_name) {
  obj <- get(obj_name)
  plot_spatial_tissue_overview(
    obj,
    tissue_name = tissue_name,
    file = file.path("summary/Fig1_overview", file_name),
    tissue_color = tissue_cols[tissue_name]
  )
})

## Fig. 1C: all-tissue UMAP overview ----
obj_list_clean <- lapply(obj_list, remove_images_for_merge)

combined_obj <- merge(
  obj_list_clean[[1]],
  y = obj_list_clean[-1],
  add.cell.ids = names(obj_list_clean),
  project = "Poplar_spatial_atlas"
)

DefaultAssay(combined_obj) <- "Spatial"
combined_obj[["Spatial"]] <- JoinLayers(combined_obj[["Spatial"]])
combined_obj <- add_organelle_qc(combined_obj, assay = "Spatial")

combined_obj <- NormalizeData(combined_obj, assay = "Spatial", verbose = FALSE)
combined_obj <- FindVariableFeatures(combined_obj, assay = "Spatial", nfeatures = 3000, verbose = FALSE)
combined_obj <- ScaleData(combined_obj, assay = "Spatial", verbose = FALSE)
combined_obj <- RunPCA(combined_obj, assay = "Spatial", npcs = 50, verbose = FALSE)
combined_obj <- RunUMAP(combined_obj, reduction = "pca", dims = 1:30)
combined_obj$tissue <- factor(combined_obj$tissue, levels = names(tissue_cols))
p_umap_tissue <- DimPlot_scCustom(combined_obj, reduction = "umap", group.by = "tissue", label = TRUE, repel = TRUE, order = TRUE, figure_plot = TRUE, colors_use = tissue_cols) +
  ggtitle("All-tissue UMAP")

ggsave("summary/Fig1_overview/Fig1C_all_tissue_UMAP.pdf", p_umap_tissue, width = 8, height = 7)

## Fig. 1C v2: all-tissue UMAP overview with Harmony integration ----
library(harmony)

obj_list_clean <- lapply(obj_list, remove_images_for_merge)
combined_obj_harmony <- merge(obj_list_clean[[1]], y = obj_list_clean[-1], add.cell.ids = names(obj_list_clean), project = "Poplar_spatial_atlas")

DefaultAssay(combined_obj_harmony) <- "Spatial"
combined_obj_harmony[["Spatial"]] <- JoinLayers(combined_obj_harmony[["Spatial"]])
combined_obj_harmony <- add_organelle_qc(combined_obj_harmony, assay = "Spatial")

combined_obj_harmony$batch <- if ("parent_library" %in% colnames(combined_obj_harmony@meta.data)) {
  combined_obj_harmony$parent_library
} else {
  combined_obj_harmony$orig.ident
}
combined_obj_harmony$batch[is.na(combined_obj_harmony$batch) | combined_obj_harmony$batch == ""] <-
  combined_obj_harmony$orig.ident[is.na(combined_obj_harmony$batch) | combined_obj_harmony$batch == ""]

combined_obj_harmony <- NormalizeData(combined_obj_harmony, assay = "Spatial", verbose = FALSE)
combined_obj_harmony <- FindVariableFeatures(combined_obj_harmony, assay = "Spatial", nfeatures = 3000, verbose = FALSE)
combined_obj_harmony <- ScaleData(combined_obj_harmony, assay = "Spatial", verbose = FALSE)
combined_obj_harmony <- RunPCA(combined_obj_harmony, assay = "Spatial", npcs = 50, verbose = FALSE)

fig1c_harmony_theta <- 3
fig1c_harmony_lambda <- 1
fig1c_harmony_dims <- 1:30

combined_obj_harmony <- RunHarmony(
  combined_obj_harmony,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use = fig1c_harmony_dims,
  theta = fig1c_harmony_theta,
  lambda = fig1c_harmony_lambda
)

fig1c_umap_neighbors <- 50
fig1c_umap_min_dist <- 0.4

combined_obj_harmony <- RunUMAP(
  combined_obj_harmony,
  reduction = "harmony",
  dims = fig1c_harmony_dims,
  n.neighbors = fig1c_umap_neighbors,
  min.dist = fig1c_umap_min_dist
)
combined_obj_harmony$tissue <- factor(combined_obj_harmony$tissue, levels = names(tissue_cols))
p_umap_tissue_harmony <- DimPlot_scCustom(combined_obj_harmony, reduction = "umap", group.by = "tissue", label = TRUE, repel = TRUE, order = TRUE, figure_plot = TRUE, colors_use = tissue_cols) +
  ggtitle("All-tissue UMAP, Harmony-integrated")

ggsave("summary/Fig1_overview/Fig1C_all_tissue_UMAP_harmony_v2.1.pdf", p_umap_tissue_harmony, width = 8, height = 7)

combined_obj_harmony$tissue_merged_petiole <- as.character(combined_obj_harmony$tissue)
combined_obj_harmony$tissue_merged_petiole[combined_obj_harmony$tissue_merged_petiole %in% c("Petiole cross", "Petiole longitudinal")] <- "Petiole"
combined_obj_harmony$tissue_merged_petiole <- factor(combined_obj_harmony$tissue_merged_petiole, levels = c("Axillary bud", "Petiole", "SAM", "Stem"))
p_umap_tissue_harmony_merged_petiole <- DimPlot_scCustom(combined_obj_harmony, reduction = "umap", group.by = "tissue_merged_petiole", label = TRUE, repel = TRUE, order = TRUE, figure_plot = TRUE, colors_use = tissue_cols_merged_petiole) +
  ggtitle("All-tissue UMAP, Harmony-integrated")

ggsave("summary/Fig1_overview/Fig1C_all_tissue_UMAP_harmony_merged_petiole.pdf", p_umap_tissue_harmony_merged_petiole, width = 8, height = 7)

# Fig. 1D: QC across tissues
p_qc1 <- VlnPlot_scCustom(combined_obj, features = "nCount_Spatial", group.by = "tissue", plot_boxplot = TRUE) + NoLegend()
p_qc2 <- VlnPlot_scCustom(combined_obj, features = "nFeature_Spatial", group.by = "tissue", plot_boxplot = TRUE) + NoLegend()
save_pdf(p_qc1 + p_qc2, "summary/Fig1_overview/Fig1D_QC_all_tissues.pdf", width = 13, height = 14)

## Fig. 1D v2: QC across tissues, log-scaled ----
tissue_order <- c("SAM", "Axillary bud", "Stem", "Petiole cross", "Petiole longitudinal")
qc_complexity_df <- combined_obj@meta.data %>%
  dplyr::select(tissue, nCount_Spatial, nFeature_Spatial) %>%
  tidyr::pivot_longer(c(nCount_Spatial, nFeature_Spatial), names_to = "metric", values_to = "value") %>%
  mutate(
    tissue = factor(tissue, levels = tissue_order),
    metric = recode(metric, nFeature_Spatial = "Detected genes per spot", nCount_Spatial = "UMI counts per spot"),
    metric = factor(metric, levels = c("Detected genes per spot", "UMI counts per spot"))
  )

qc_organelle_df <- combined_obj@meta.data %>%
  dplyr::select(tissue, percent_organelle) %>%
  tidyr::pivot_longer(percent_organelle, names_to = "metric", values_to = "value") %>%
  mutate(
    tissue = factor(tissue, levels = tissue_order),
    metric = recode(
      metric,
      percent_organelle = "ChrPt + ChrMt UMIs without rrnL/rrnS (%)"
    ),
    metric = factor(metric, levels = "ChrPt + ChrMt UMIs without rrnL/rrnS (%)")
  )

qc_organelle_summary <- combined_obj@meta.data %>%
  group_by(tissue) %>%
  summarise(
    n_spots = dplyr::n(),
    median_chrpt_percent = median(percent_chrpt, na.rm = TRUE),
    median_chrmt_all_percent = median(percent_chrmt_all, na.rm = TRUE),
    median_chrmt_rrnLS_percent = median(percent_chrmt_rrnLS, na.rm = TRUE),
    median_chrmt_percent = median(percent_chrmt, na.rm = TRUE),
    median_organelle_all_percent = median(percent_organelle_all, na.rm = TRUE),
    median_organelle_percent = median(percent_organelle, na.rm = TRUE),
    mean_chrpt_percent = mean(percent_chrpt, na.rm = TRUE),
    mean_chrmt_all_percent = mean(percent_chrmt_all, na.rm = TRUE),
    mean_chrmt_rrnLS_percent = mean(percent_chrmt_rrnLS, na.rm = TRUE),
    mean_chrmt_percent = mean(percent_chrmt, na.rm = TRUE),
    mean_organelle_all_percent = mean(percent_organelle_all, na.rm = TRUE),
    mean_organelle_percent = mean(percent_organelle, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(tissue = factor(tissue, levels = tissue_order)) %>%
  arrange(tissue)

write.csv(qc_organelle_summary, "summary/Fig1_overview/Fig1D_organelle_percent_summary.csv", row.names = FALSE)
write.csv(
  data.frame(
    excluded_feature = c("gene-rrnL", "gene-rrnS"),
    gene_name = c("rrnL", "rrnS"),
    metric = "percent_chrmt and percent_organelle",
    reason = "Dominant mitochondrial rRNA features excluded from plotted organelle QC metric"
  ),
  "summary/Fig1_overview/Fig1D_organelle_rrnLS_excluded_features.csv",
  row.names = FALSE
)

p_qc_v2 <- ggplot(qc_complexity_df, aes(x = tissue, y = value, fill = tissue)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "black") +
  geom_boxplot(width = 0.12, outlier.size = 0.12, linewidth = 0.25, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_cols) +
  scale_y_log10() +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(y = "Value per spot, log10 scale")

p_qc_v2_horizontal <- ggplot(qc_complexity_df, aes(x = value, y = tissue, fill = tissue)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "black") +
  geom_boxplot(width = 0.12, outlier.size = 0.12, linewidth = 0.25, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_cols) +
  scale_x_log10() +
  facet_wrap(~metric, ncol = 1, scales = "free_x") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(x = "Value per spot, log10 scale")

p_qc_organelle <- ggplot(qc_organelle_df, aes(x = tissue, y = value, fill = tissue)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "black") +
  geom_boxplot(width = 0.12, outlier.size = 0.12, linewidth = 0.25, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_cols) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(y = "UMIs per spot (%)")

p_qc_organelle_horizontal <- ggplot(qc_organelle_df, aes(x = value, y = tissue, fill = tissue)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "black") +
  geom_boxplot(width = 0.12, outlier.size = 0.12, linewidth = 0.25, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_cols) +
  facet_wrap(~metric, ncol = 1, scales = "free_x") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(x = "UMIs per spot (%)")

p_qc_gene_main <- ggplot(combined_obj@meta.data, aes(x = factor(tissue, levels = tissue_order), y = nFeature_Spatial, fill = tissue)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "black") +
  geom_boxplot(width = 0.12, outlier.size = 0.12, linewidth = 0.25, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_cols) +
  scale_y_log10() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
  ) +
  labs(title = "Detected genes per spot", y = "Gene number, log10 scale")

p_qc_umi_main <- ggplot(combined_obj@meta.data, aes(x = factor(tissue, levels = tissue_order), y = nCount_Spatial, fill = tissue)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "black") +
  geom_boxplot(width = 0.12, outlier.size = 0.12, linewidth = 0.25, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_cols) +
  scale_y_log10() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
  ) +
  labs(title = "UMI counts per spot", y = "UMI number, log10 scale")

p_qc_organelle_main <- ggplot(combined_obj@meta.data, aes(x = factor(tissue, levels = tissue_order), y = percent_organelle, fill = tissue)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "black") +
  geom_boxplot(width = 0.12, outlier.size = 0.12, linewidth = 0.25, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_cols) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
  ) +
  labs(title = "ChrPt + ChrMt UMIs without rrnL/rrnS (%)", y = "mt + cp UMIs per spot (%)")

p_qc_with_organelle <- p_qc_gene_main / p_qc_umi_main / p_qc_organelle_main + plot_layout(heights = c(1, 1, 1))

qc_tissue_level_df <- combined_obj@meta.data %>%
  mutate(
    tissue_level = as.character(tissue),
    tissue_level = ifelse(tissue_level %in% c("Petiole cross", "Petiole longitudinal"), "Petiole", tissue_level),
    tissue_level = factor(tissue_level, levels = c("SAM", "Axillary bud", "Stem", "Petiole"))
  )

qc_metrics_section_table <- combined_obj@meta.data %>%
  mutate(grouping_level = "section", tissue_label = as.character(tissue)) %>%
  group_by(grouping_level, tissue_label) %>%
  summarise(
    n_spots = dplyr::n(),
    total_umis = sum(nCount_Spatial, na.rm = TRUE),
    median_detected_genes = median(nFeature_Spatial, na.rm = TRUE),
    q1_detected_genes = quantile(nFeature_Spatial, 0.25, na.rm = TRUE),
    q3_detected_genes = quantile(nFeature_Spatial, 0.75, na.rm = TRUE),
    mean_detected_genes = mean(nFeature_Spatial, na.rm = TRUE),
    median_umi_counts = median(nCount_Spatial, na.rm = TRUE),
    q1_umi_counts = quantile(nCount_Spatial, 0.25, na.rm = TRUE),
    q3_umi_counts = quantile(nCount_Spatial, 0.75, na.rm = TRUE),
    mean_umi_counts = mean(nCount_Spatial, na.rm = TRUE),
    median_chrpt_percent = median(percent_chrpt, na.rm = TRUE),
    median_chrmt_percent_without_rrnLS = median(percent_chrmt, na.rm = TRUE),
    median_chrpt_chrmt_percent_without_rrnLS = median(percent_organelle, na.rm = TRUE),
    median_chrmt_percent_all = median(percent_chrmt_all, na.rm = TRUE),
    median_chrmt_rrnLS_percent = median(percent_chrmt_rrnLS, na.rm = TRUE),
    median_chrpt_chrmt_percent_all = median(percent_organelle_all, na.rm = TRUE),
    mean_chrpt_chrmt_percent_without_rrnLS = mean(percent_organelle, na.rm = TRUE),
    .groups = "drop"
  )

qc_metrics_tissue_table <- qc_tissue_level_df %>%
  mutate(grouping_level = "tissue", tissue_label = as.character(tissue_level)) %>%
  group_by(grouping_level, tissue_label) %>%
  summarise(
    n_spots = dplyr::n(),
    total_umis = sum(nCount_Spatial, na.rm = TRUE),
    median_detected_genes = median(nFeature_Spatial, na.rm = TRUE),
    q1_detected_genes = quantile(nFeature_Spatial, 0.25, na.rm = TRUE),
    q3_detected_genes = quantile(nFeature_Spatial, 0.75, na.rm = TRUE),
    mean_detected_genes = mean(nFeature_Spatial, na.rm = TRUE),
    median_umi_counts = median(nCount_Spatial, na.rm = TRUE),
    q1_umi_counts = quantile(nCount_Spatial, 0.25, na.rm = TRUE),
    q3_umi_counts = quantile(nCount_Spatial, 0.75, na.rm = TRUE),
    mean_umi_counts = mean(nCount_Spatial, na.rm = TRUE),
    median_chrpt_percent = median(percent_chrpt, na.rm = TRUE),
    median_chrmt_percent_without_rrnLS = median(percent_chrmt, na.rm = TRUE),
    median_chrpt_chrmt_percent_without_rrnLS = median(percent_organelle, na.rm = TRUE),
    median_chrmt_percent_all = median(percent_chrmt_all, na.rm = TRUE),
    median_chrmt_rrnLS_percent = median(percent_chrmt_rrnLS, na.rm = TRUE),
    median_chrpt_chrmt_percent_all = median(percent_organelle_all, na.rm = TRUE),
    mean_chrpt_chrmt_percent_without_rrnLS = mean(percent_organelle, na.rm = TRUE),
    .groups = "drop"
  )

qc_metrics_supp_table <- bind_rows(qc_metrics_section_table, qc_metrics_tissue_table) %>%
  mutate(
    grouping_level = factor(grouping_level, levels = c("section", "tissue")),
    tissue_label = factor(tissue_label, levels = c(tissue_order, "Petiole"))
  ) %>%
  arrange(grouping_level, tissue_label)

write.csv(qc_metrics_supp_table, "summary/Supplementary/Supp_Table_QC_metrics_by_tissue.csv", row.names = FALSE)

p_qc_gene_tissue_level <- ggplot(qc_tissue_level_df, aes(x = tissue_level, y = nFeature_Spatial, fill = tissue_level)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "black") +
  geom_boxplot(width = 0.12, outlier.size = 0.12, linewidth = 0.25, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_cols_merged_petiole) +
  scale_y_log10() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
  ) +
  labs(title = "Detected genes per spot", y = "Gene number, log10 scale")

p_qc_umi_tissue_level <- ggplot(qc_tissue_level_df, aes(x = tissue_level, y = nCount_Spatial, fill = tissue_level)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "black") +
  geom_boxplot(width = 0.12, outlier.size = 0.12, linewidth = 0.25, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_cols_merged_petiole) +
  scale_y_log10() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
  ) +
  labs(title = "UMI counts per spot", y = "UMI number, log10 scale")

p_qc_organelle_tissue_level <- ggplot(qc_tissue_level_df, aes(x = tissue_level, y = percent_organelle, fill = tissue_level)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "black") +
  geom_boxplot(width = 0.12, outlier.size = 0.12, linewidth = 0.25, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_cols_merged_petiole) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
  ) +
  labs(title = "ChrPt + ChrMt UMIs without rrnL/rrnS (%)", y = "mt + cp UMIs per spot (%)")

p_qc_with_organelle_tissue_level <- p_qc_gene_tissue_level / p_qc_umi_tissue_level / p_qc_organelle_tissue_level + plot_layout(heights = c(1, 1, 1))

ggsave("summary/Fig1_overview/Fig1D_QC_all_tissues_v2_logscale.pdf", p_qc_v2, width = 8, height = 7)
ggsave("summary/Fig1_overview/Fig1D_QC_all_tissues_v2_horizontal_logscale.pdf", p_qc_v2_horizontal, width = 6, height = 8)
ggsave("summary/Fig1_overview/Fig1D_main_QC_gene_UMI_organelle_rrnLS_excluded.pdf", p_qc_with_organelle, width = 8, height = 9)
ggsave("summary/Fig1_overview/Fig1D_organelle_percent_by_tissue.pdf", p_qc_organelle, width = 8, height = 7)
ggsave("summary/Fig1_overview/Fig1D_organelle_percent_by_tissue_horizontal.pdf", p_qc_organelle_horizontal, width = 6, height = 8)
ggsave("summary/Fig1_overview/Fig1D_QC_all_tissues_with_organelle_percent.pdf", p_qc_with_organelle, width = 8, height = 9)
ggsave("summary/Fig1_overview/Fig1D_QC_tissue_level_with_organelle_percent.pdf", p_qc_with_organelle_tissue_level, width = 8, height = 9)

if (identical(Sys.getenv("RUN_FIG1_QC_ONLY"), "1")) {
  message("RUN_FIG1_QC_ONLY complete after Figure 1 QC outputs.")
  quit(save = "no", status = 0)
}

# --------------------------------------- #
# 2. Cell/domain annotation figure ####
# --------------------------------------- #

#### Update SAM annotation for plotting 
# [2026.7.24 REVISION] Re-mapped to the lower-leaf-cleaned SAM object
#   (sobj_sam_split_cleaned_res0.4_annotated_v2_2026.7.qs, 5,267 spots, res 0.4 with the finer
#   5_0/5_1 and 8_0/8_1/8_2 sub-cluster splits). Labels are the user's manual
#   "marker eye-test" calls after cleaning: vasculature is now resolved into
#   Xylem (3) vs Phloem (4); cluster 10 = Mesophyll; 8_0/9 = Proliferating;
#   8_2 = Shoot meristem. Old (pre-cleaning) 0-10 flat mapping kept below in comments.
sam_anno_v2 <- tibble::tribble(
  ~cluster_label, ~celltype_predicted, ~predicted_score, ~celltype_final, ~celltypes_v2,
  "0",   "Pith",        0.44, "Pith",            "Pith",
  "1",   "Cortex",      0.28, "Cortex",          "Cortex",
  "2",   "Epidermis",   0.19, "Epidermis",       "Epidermis",
  "3",   "Vasculature", 0.19, "Xylem",           "Xylem",
  "4",   "Vasculature", 0.11, "Phloem",          "Phloem",
  "5_0", "Epidermis",   0.16, "Epidermis",       "Epidermis",
  "5_1", "Cortex",      0.16, "Cortex",          "Cortex",
  "6",   "Epidermis",   0.16, "Leaf primordium", "Leaf primordium",
  "7",   "Epidermis",   0.16, "Epidermis",       "Epidermis",
  "8_0", "Proliferating", 0.11, "Proliferating", "Proliferating",
  "8_1", "Epidermis",   0.11, "Leaf primordium", "Leaf primordium",
  "8_2", "Meristem",    0.05, "Shoot meristem",  "Shoot meristem",
  "9",   "Proliferating", 0.05, "Proliferating", "Proliferating",
  "10",  "Mesophyll",   0.09, "Mesophyll",       "Mesophyll"
)
# --- OLD pre-cleaning mapping (superseded 2026.7.24) ---
#   "0"=Pith, "1"=Vasculature, "2"=Cortex, "3/5/6_1/7/8"=Leaf primordium1-5,
#   "4/10"=Epidermis, "6_0"=Axillary meristem, "9"=Meristem

## keep original cluster number and v1 annotation
sobj_sam_ann$sam_cluster_label <- as.character(sobj_sam_ann$seurat_clusters)
sobj_sam_ann$celltypes_v1 <- sobj_sam_ann$celltypes

## add v2 annotation
sam_v2_map <- setNames(sam_anno_v2$celltypes_v2, sam_anno_v2$cluster_label)
sobj_sam_ann$celltypes_v2 <- unname(sam_v2_map[sobj_sam_ann$sam_cluster_label])
sobj_sam_ann$celltype <- sobj_sam_ann$celltypes_v2
Idents(sobj_sam_ann) <- "celltype"
qsave(sobj_sam_ann, "saved_obj/sobj_sam_split_cleaned_res0.4_annotated_v2_2026.7.qs")

###### Update Bud annotation for plotting
bud_anno_v2 <- tibble::tribble(
  ~cluster_label, ~celltype_predicted, ~predicted_score, ~celltype_final, ~celltypes_v2,
  "0",  "Cortex",        0.329545653, "Cortex1",        "Cortex",
  "1",  "Cortex",        0.690734351, "Cortex2",        "Cortex",
  "2",  "Cortex",        0.470095123, "Epidermis",      "Cortex",
  "3",  "Sieve Element", 0.271646914, "Vasculature1",   "Vasculature",
  "4",  "Cortex",        0.351215040, "Vasculature2",   "Phloem",
  "5",  "Cortex",        0.592866054, "Epidermis",      "Epidermis",
  "6",  "Cortex",        0.473046651, "Vasculature3",   "Xylem",
  "7",  "Cortex",        0.493810934, "Cortex3",        "Cortex",
  "8",  "Epidermis",     0.250341892, "Bud primordium", "Bud primordium",
  "9",  "Epidermis",     0.687767072, "Bud scale",      "Bud scale",
  "10", "Cortex",        0.585253275, "Cortex4",        "Cortex",
  "11", "Pith",          0.144730942, "Vasculature4",   "Bud procambium",
  "12", "Epidermis",     0.289340025, "Epidermis",      "Epidermis",
  "13", "Cortex",        0.581586710, "Cortex5",        "Cortex",
  "14", "Cortex",        0.458520051, "Epidermis",      "Epidermis"
)
## keep original cluster number and v1 annotation
sobj_bud_ann$bud_cluster_label <- as.character(sobj_bud_ann$seurat_clusters)
sobj_bud_ann$celltypes_v1 <- sobj_bud_ann$celltypes

## add v2 annotation
bud_v2_map <- setNames(bud_anno_v2$celltypes_v2, bud_anno_v2$cluster_label)
sobj_bud_ann$celltypes_v2 <- unname(bud_v2_map[sobj_bud_ann$bud_cluster_label])

## fallback: keep v1 label if any cluster is missing from the table
sobj_bud_ann$celltypes_v2[is.na(sobj_bud_ann$celltypes_v2)] <- sobj_bud_ann$celltypes_v1[is.na(sobj_bud_ann$celltypes_v2)]

table(sobj_bud_ann$bud_cluster_label, sobj_bud_ann$celltypes_v1)
table(sobj_bud_ann$celltypes_v2)
sobj_bud_ann$celltype <- sobj_bud_ann$celltypes_v2
Idents(sobj_bud_ann) <- "celltype"
qsave(sobj_bud_ann, "saved_obj/sobj_bud_res0.5_final_v2_2026.7.24.qs")
sobj_stem_cross_ann$celltype <- as.character(sobj_stem_cross_ann$celltypes)
Idents(sobj_stem_cross_ann) <- "celltype"
qsave(sobj_stem_cross_ann, "saved_obj/sobj_stem_cross_AB_res0.6_final_v2_Seurat_RCTD_2026.7.26.qs")
sobj_petiole_cross_ann$celltype <- as.character(sobj_petiole_cross_ann$celltypes)
Idents(sobj_petiole_cross_ann) <- "celltype"
qsave(sobj_petiole_cross_ann, "saved_obj/sobj_petiole_cross_res0.35_final_v1_2026.7.24.qs")
sobj_petiole_long_ann$celltype <- as.character(sobj_petiole_long_ann$celltypes)
Idents(sobj_petiole_long_ann) <- "celltype"
qsave(sobj_petiole_long_ann, "saved_obj/sobj_petiole_longitudinal_res0.15_final_v1_2026.7.24.qs")

final_celltype_audit <- bind_rows(
  tibble(tissue = "SAM", celltype = as.character(sobj_sam_ann$celltype)),
  tibble(tissue = "Axillary bud", celltype = as.character(sobj_bud_ann$celltype)),
  tibble(tissue = "Stem", celltype = as.character(sobj_stem_cross_ann$celltype)),
  tibble(tissue = "Petiole cross", celltype = as.character(sobj_petiole_cross_ann$celltype)),
  tibble(tissue = "Petiole longitudinal", celltype = as.character(sobj_petiole_long_ann$celltype))
) %>% count(tissue, celltype, name = "n_spots") %>% arrange(tissue, celltype)
write_csv(final_celltype_audit, "summary/Fig2_annotation/final_celltype_label_audit.csv")

# Fig. 2A: annotated spatial maps for each tissue
plot_spatial_annotation(sobj_bud_ann, "Axillary bud annotation", "summary/Fig2_annotation/Fig2A_bud_spatial_annotation.pdf", group_col = "celltype", width = 16, height = 9)
plot_spatial_annotation(sobj_stem_cross_ann, "Stem annotation", "summary/Fig2_annotation/Fig2A_stem_spatial_annotation.pdf", group_col = "celltype", width = 16, height = 9)
plot_spatial_annotation(sobj_sam_ann, "SAM annotation", "summary/Fig2_annotation/Fig2A_sam_spatial_annotation_v2.pdf", group_col = "celltype", width = 16, height = 9)
plot_spatial_annotation(sobj_petiole_cross_ann, "Petiole cross annotation", "summary/Fig2_annotation/Fig2A_petiole_cross_spatial_annotation.pdf", group_col = "celltype", width = 16, height = 9)
plot_spatial_annotation(sobj_petiole_long_ann, "Petiole longitudinal annotation", "summary/Fig2_annotation/Fig2A_petiole_long_spatial_annotation.pdf", group_col = "celltype", width = 16, height = 9)

# Fig. 2B: UMAP annotation for each tissue
plot_umap_annotation(sobj_bud_ann, "Axillary bud UMAP", "summary/Fig2_annotation/Fig2B_bud_UMAP_annotation.pdf", group_col = "celltype")
plot_umap_annotation(sobj_stem_cross_ann, "Stem UMAP", "summary/Fig2_annotation/Fig2B_stem_UMAP_annotation.pdf", group_col = "celltype")
plot_umap_annotation(sobj_sam_ann, "SAM UMAP", "summary/Fig2_annotation/Fig2B_sam_UMAP_annotation.pdf", group_col = "celltype")
plot_umap_annotation(sobj_petiole_cross_ann, "Petiole cross UMAP", "summary/Fig2_annotation/Fig2B_petiole_cross_UMAP_annotation.pdf", group_col = "celltype")
plot_umap_annotation(sobj_petiole_long_ann, "Petiole longitudinal UMAP", "summary/Fig2_annotation/Fig2B_petiole_long_UMAP_annotation.pdf", group_col = "celltype")


## Plot marker heatmap again for updated object:
# Heatmap for cluster interpretation
if (exists("make_marker_heatmap") && exists("mkr_list_use")) {
  make_marker_heatmap(
    sobj = sobj_sam_ann,
    marker_df = mkr_list_use,
    out_pdf = "marker/marker_heatmap_split_celltype_SAM_annotation_V2.pdf"
  )

  make_marker_heatmap(
    sobj = sobj_bud_ann,
    marker_df = mkr_list_use,
    out_pdf = "marker/marker_heatmap_split_celltype_BUD_annotation_V2.pdf"
  )

  make_marker_heatmap(
    sobj = sobj_petiole_cross_ann,
    marker_df = mkr_list_use,
    out_pdf = "marker/marker_heatmap_split_celltype_PETIOLE_CROSS_annotation_V2.pdf"
  )
} else {
  message("Optional old marker heatmaps skipped: make_marker_heatmap() or mkr_list_use is not available.")
}




# Fig. 2C: marker dotplots
# These marker files are generated in tissue-specific pipelines. Update names if needed.
bud_markers <- safe_read_csv("marker/BUD_all_markers_by_celltype_PtXaOnly.csv")
stem_markers <- safe_read_csv("marker/STEM_long_CD_all_markers_by_celltype_PtXaOnly.csv")
sam_markers <- safe_read_csv("marker/SAM_all_markers_by_celltype.csv")
petiole_cross_markers <- safe_read_csv("marker/PETIOLE_CROSS_all_markers_by_celltype_PtXaOnly.csv")
petiole_long_markers <- safe_read_csv("marker/PETIOLE_LONG_all_markers_by_celltype_PtXaOnly.csv")

plot_marker_dotplot(sobj_bud_ann, bud_markers, "summary/Fig2_annotation/Fig2C_bud_top_marker_dotplot.pdf")
plot_marker_dotplot(sobj_stem_cross_ann, stem_markers, "summary/Fig2_annotation/Fig2C_stem_top_marker_dotplot.pdf")
plot_marker_dotplot(sobj_sam_ann, sam_markers, "summary/Fig2_annotation/Fig2C_sam_top_marker_dotplot.pdf")
plot_marker_dotplot(sobj_petiole_cross_ann, petiole_cross_markers, "summary/Fig2_annotation/Fig2C_petiole_cross_top_marker_dotplot.pdf")
plot_marker_dotplot(sobj_petiole_long_ann, petiole_long_markers, "summary/Fig2_annotation/Fig2C_petiole_long_top_marker_dotplot.pdf")



# Fig. 2C: reviewed curated marker dotplots for selected representative tissues ----
# The source marker sheet is used as an audit table; the plotted marker set below
# keeps two high-confidence markers per annotated tissue domain wherever possible.

marker_source_file <- c(
  "updated_markerlist_poplar_2.8.26.csv",
  "marker/updated_markerlist_poplar_2.8.26.csv",
  "tables/updated_markerlist_poplar_2.8.26.csv"
) %>%
  purrr::keep(file.exists) %>%
  purrr::pluck(1, .default = NA_character_)

updated_marker_source <- if (!is.na(marker_source_file)) {
  safe_read_csv(marker_source_file)
} else {
  message("Missing updated_markerlist_poplar_2.8.26.csv; using reviewed marker table only.")
  NULL
}

fig2c_reviewed_markers <- tibble::tribble(
  ~tissue, ~celltype, ~marker_name, ~gene_id, ~marker_role, ~plot_set, ~notes,
  # [2026.7.24 REVISION] SAM markers re-curated against the lower-leaf-cleaned
  #   object + a per-marker specificity audit (z-scored avg expr across the 9
  #   celltypes_v2 domains; "peaks in intended domain?" test). plot_set = "main"
  #   are canonical genes that DO peak in their domain; "supplement" documents
  #   classic-but-non-specific or redundant markers (kept traceable, off-plot).
  #   CAVEAT: Cortex has no clean positive marker — cortex & pith are both ground
  #   parenchyma and overlap transcriptionally; PMEAMT/PSBXb peak in Pith. Marked
  #   main for completeness but flagged in notes.
  "SAM", "Shoot meristem", "KNAT2", "PtXaAlbH.10G030700.v5.1", "meristem identity", "main", "Class-I KNOX homeobox; canonical SAM/indeterminacy marker. z=2.27 in meristem.",
  "SAM", "Shoot meristem", "BLH8", "PtXaAlbH.04G165800.v5.1", "meristem identity", "main", "BELL1-like homeobox (PNY); meristem/organ-boundary. z=1.85.",
  "SAM", "Shoot meristem", "RPL", "PtXaAlbH.10G156800.v5.1", "meristem identity", "main", "REPLUMLESS/BLH9 homeobox; strongest meristem marker. z=2.65.",
  "SAM", "Shoot meristem", "H1", "PtXaAlbH.18G040700.v5.1", "proliferation state", "supplement", "Linker histone; meristem-enriched but generic.",
  "SAM", "Shoot meristem", "MCM7", "PtXaTreH.14G094300.v5.1", "proliferation state", "supplement", "DNA-replication licensing; overlaps proliferating.",
  "SAM", "Proliferating", "CYCA2;3", "PtXaAlbH.01G152100.v5.1", "cell cycle", "main", "A2 cyclin; canonical G2-M cell-cycle marker. z=2.01.",
  "SAM", "Proliferating", "PCNA2", "PtXaTreH.01G204900.v5.1", "cell cycle", "main", "Proliferating cell nuclear antigen; S-phase. z=1.88.",
  "SAM", "Proliferating", "GAMMA-H2AX", "PtXaTreH.05G027800.v5.1", "cell cycle", "main", "Replication-linked histone variant. z=2.23.",
  "SAM", "Proliferating", "GL1", "PtXaTreH.12G058600.v5.1", "trichome initiation", "supplement", "GLABRA1; kept for trichome story (Fig4/5), not domain ID.",
  "SAM", "Leaf primordium", "HIC", "PtXaTreH.02G151000.v5.1", "primordium emergence", "main", "Epidermal/guard-cell; strong LP-emergence marker (71% spots, z=2.05).",
  "SAM", "Leaf primordium", "LTL1", "PtXaAlbH.19G028400.v5.1", "epidermis/cuticle", "main", "Lipid-transfer-like; LP epidermal. z=2.07.",
  "SAM", "Leaf primordium", "TRIPTYCHON", "PtXaTreH.15G015900.v5.1", "trichome patterning", "supplement", "Peaks in proliferating, not LP; documented only.",
  "SAM", "Epidermis", "PDF1", "PtXaTreH.02G052300.v5.1", "protoderm identity", "main", "PROTODERMAL FACTOR1; canonical protoderm/epidermis. z=1.98.",
  "SAM", "Epidermis", "KCS2", "PtXaAlbH.10G060600.v5.1", "epidermis/cuticle", "main", "3-ketoacyl-CoA synthase; cuticular wax. z=1.85.",
  "SAM", "Epidermis", "DCR", "PtXaAlbH.15G116500.v5.1", "cutin biosynthesis", "main", "DEFECTIVE IN CUTICULAR RIDGES; cutin. z=1.63.",
  "SAM", "Epidermis", "GDSL-lipase", "PtXaTreH.04G052500.v5.1", "epidermis/cuticle", "supplement", "Epidermal but shared with LP; redundant.",
  "SAM", "Epidermis", "LPTG1", "PtXaTreH.01G047400.v5.1", "epidermis/cuticle", "supplement", "GPI-anchored LTP; redundant with PDF1/KCS2.",
  "SAM", "Epidermis", "RWP1-like", "PtXaAlbH.13G067400.v5.1", "epidermis/cuticle", "supplement", "Fails specificity (peaks in pith).",
  "SAM", "Epidermis", "VEP1", "PtXaAlbH.02G104000.v5.1", "epidermis/cuticle", "supplement", "Fails specificity (peaks in pith).",
  "SAM", "Mesophyll", "CRR23", "PtXaAlbH.10G085700.v5.1", "photosynthetic identity", "main", "Chlororespiratory NDH; photosynthetic mesophyll. z=2.16.",
  "SAM", "Mesophyll", "CER5", "PtXaTreH.01G212400.v5.1", "wax export", "main", "ECERIFERUM5/ABCG12 wax export; mesophyll-enriched here. z=2.32.",
  "SAM", "Mesophyll", "CA1", "PtXaTreH.03G091500.v5.1", "photosynthetic identity", "supplement", "Classic mesophyll carbonic anhydrase; peaks proliferating here.",
  "SAM", "Mesophyll", "FIL", "PtXaTreH.02G123100.v5.1", "abaxial identity", "supplement", "FILAMENTOUS FLOWER; peaks proliferating here.",
  "SAM", "Cortex", "PMEAMT", "PtXaAlbH.15G031100.v5.1", "domain-enriched", "main", "Pectin methylesterase; best-available cortex (CAVEAT: peaks pith, z_cortex=0.69).",
  "SAM", "Cortex", "PSBXb", "PtXaTreH.06G124900.v5.1", "domain-enriched", "main", "Broad ground tissue; 92% cortex spots (CAVEAT: peaks pith).",
  "SAM", "Cortex", "PNP", "PtXaAlbH.13G084300.v5.1", "domain-enriched", "supplement", "Plant natriuretic peptide; weak/pith-shared.",
  "SAM", "Pith", "WRKY12", "PtXaTreH.14G037200.v5.1", "identity", "main", "Canonical pith secondary-wall repressor. z=1.92.",
  "SAM", "Pith", "ATHB13", "PtXaTreH.10G076500.v5.1", "identity", "main", "HD-Zip I; pith-enriched. z=2.28.",
  "SAM", "Pith", "WRKY13", "PtXaTreH.07G066400.v5.1", "identity", "supplement", "WRKY13; pith but low detection (7% spots).",
  "SAM", "Xylem", "NST1", "PtXaTreH.14G081200.v5.1", "secondary wall master TF", "main", "NAC secondary-wall master TF; canonical xylem/fiber. z=2.63.",
  "SAM", "Xylem", "XCP1", "PtXaAlbH.04G160300.v5.1", "vessel differentiation", "main", "Xylem cysteine peptidase; tracheary-element PCD. z=2.65.",
  "SAM", "Xylem", "FLA12", "PtXaTreH.12G095800.v5.1", "secondary wall", "main", "Fasciclin-like AGP; secondary-wall/xylem. z=2.67.",
  "SAM", "Xylem", "MAN6", "PtXaAlbH.16G114900.v5.1", "secondary wall", "supplement", "Endo-mannanase; redundant with NST1/XCP1.",
  "SAM", "Phloem", "SEOR1", "PtXaAlbH.01G281100.v5.1", "sieve element", "main", "Sieve element occlusion-related; canonical SE. z=2.53.",
  "SAM", "Phloem", "CLE41", "PtXaTreH.02G198800.v5.1", "phloem-procambium identity", "main", "TDIF/CLE41 peptide; phloem-procambium. z=2.57.",
  "SAM", "Phloem", "AN5", "PtXaAlbH.12G091700.v5.1", "companion cell", "main", "Companion-cell associated; strong phloem. z=2.61.",
  "SAM", "Phloem", "SUT4", "PtXaTreH.02G092500.v5.1", "sucrose transport", "supplement", "Sucrose transporter; lower detection.",
  "SAM", "Phloem", "URP3", "PtXaAlbH.01G198900.v5.1", "phloem support", "supplement", "Phloem-associated; redundant.",
  "SAM", "Phloem", "GSL11", "PtXaTreH.02G050500.v5.1", "callose synthase", "supplement", "Callose synthase; phloem, low detection.",

  "Axillary bud", "Cortex", "PMEAMT", "PtXaAlbH.15G031100.v5.1", "domain-enriched", "main", "Dataset-supported cortex marker.",
  "Axillary bud", "Cortex", "PNP", "PtXaAlbH.13G084300.v5.1", "domain-enriched", "main", "Added as second cortex-supporting marker if present.",
  "Axillary bud", "Vasculature", "PXY", "PtXaTreH.01G106100.v5.1", "procambium/cambium identity", "main", "Strong vascular/procambium marker.",
  "Axillary bud", "Vasculature", "FLA12", "PtXaTreH.12G095800.v5.1", "vascular support", "main", "Useful vascular-domain support marker.",
  "Axillary bud", "Phloem", "GSL11", "PtXaAlbH.02G051200.v5.1", "phloem/sieve support", "main", "Phloem-supporting marker from hand-curated list.",
  "Axillary bud", "Phloem", "UMAMIT12", "PtXaAlbH.06G071600.v5.1", "phloem/companion support", "main", "Added from reviewed phloem marker set if present.",
  "Axillary bud", "Xylem", "LAC4", "PtXaAlbH.06G084500.v5.1", "xylem/lignification", "main", "Strong xylem-supporting marker.",
  "Axillary bud", "Xylem", "4CL", "PtXaAlbH.01G031600.v5.1", "xylem/lignification", "main", "Strong lignification/secondary-wall marker.",
  "Axillary bud", "Epidermis", "CER5", "PtXaTreH.01G212400.v5.1", "epidermis/cuticle", "main", "Good epidermal/cuticle marker.",
  "Axillary bud", "Epidermis", "KCS2", "PtXaAlbH.10G060600.v5.1", "epidermis/cuticle", "main", "Good epidermal/cuticle marker.",
  "Axillary bud", "Bud primordium", "SCL28", "PtXaAlbH.01G093300.v5.1", "proliferative primordium state", "main", "Supports proliferative primordium state.",
  "Axillary bud", "Bud primordium", "CYCB1;5", "PtXaAlbH.16G025500.v5.1", "proliferation state", "main", "State marker; interpret with anatomy.",
  "Axillary bud", "Bud scale", "DCR", "PtXaAlbH.15G116500.v5.1", "cuticle/bud-scale support", "main", "Preferred bud-scale marker from hand-curated list.",
  "Axillary bud", "Bud scale", "PEC1", "PtXaAlbH.18G054000.v5.1", "epidermal/bud-scale support", "main", "Use with DCR for bud-scale support.",
  "Axillary bud", "Bud procambium", "PXY", "PtXaTreH.01G106100.v5.1", "procambium/cambium identity", "main", "Preferred if detected in this object.",
  "Axillary bud", "Bud procambium", "ANT-3", "PtXaAlbH.14G005300.v5.1", "procambium/development", "main", "Developmental marker; interpret with spatial restriction.",

  "Petiole cross", "Epidermis", "FDH", "PtXaTreH.06G179500.v5.1", "epidermis/cuticle", "main", "Strong epidermal/cuticle marker.",
  "Petiole cross", "Epidermis", "LTL1", "PtXaAlbH.19G028400.v5.1", "epidermis/cuticle", "main", "Preferred over PEC1 for main epidermis validation.",
  "Petiole cross", "Vasculature (Phloem)", "UMAMIT12", "PtXaAlbH.06G071600.v5.1", "phloem/companion support", "main", "Good phloem-supporting marker.",
  "Petiole cross", "Vasculature (Phloem)", "OPS", "PtXaAlbH.06G084100.v5.1", "phloem/procambium support", "main", "Preferred with UMAMIT12; use SUT4 if OPS is absent.",
  "Petiole cross", "Pith (Xylem)", "LAC4", "PtXaAlbH.06G084500.v5.1", "xylem/lignification", "main", "Strong xylem-supporting marker.",
  "Petiole cross", "Pith (Xylem)", "WRKY12", "PtXaTreH.14G037200.v5.1", "pith/xylem support", "main", "Retained from hand-curated list for mixed pith/xylem domain.",
  "Petiole cross", "Inner Cortex", "PNP", "PtXaAlbH.13G084300.v5.1", "domain-enriched", "main", "Dataset-supported inner cortex marker.",
  "Petiole cross", "Inner Cortex", "PMEAMT", "PtXaAlbH.15G031100.v5.1", "domain-enriched", "main", "Dataset-supported inner cortex marker.",
  "Petiole cross", "Cortex", "PSBXb", "PtXaTreH.06G124900.v5.1", "domain-enriched", "main", "Use as cortex-domain marker only if spatially specific.",
  "Petiole cross", "Cortex", "HYDROLASE", "PtXaAlbH.18G083900.v5.1", "domain-enriched", "main", "Use as cortex-domain marker only if spatially specific.",

  "Petiole longitudinal", "Cortex1", "CLE41a", "PtXaAlbH.02G198600.v5.1", "cambium/procambium identity", "main", "Validated longitudinal-petiole Cortex1 marker.",
  "Petiole longitudinal", "Cortex1", "CYP704A1", "PtXaTreH.14G053100.v5.1", "lipid metabolism", "main", "Validated longitudinal-petiole Cortex1 marker.",
  "Petiole longitudinal", "Cortex2", "CAB3", "PtXaTreH.05G186800.v5.1", "photosynthetic tissue", "main", "Validated longitudinal-petiole Cortex2 marker.",
  "Petiole longitudinal", "Epidermis", "LTL1", "PtXaTreH.19G028800.v5.1", "epidermis/cuticle", "main", "Validated longitudinal-petiole epidermis marker.",
  "Petiole longitudinal", "Epidermis", "LTPG2", "PtXaTreH.09G132100.v5.1", "epidermis/cuticle", "main", "Validated longitudinal-petiole epidermis marker.",
  "Petiole longitudinal", "Phloem", "PP2-A10", "PtXaAlbH.05G153700.v5.1", "phloem", "main", "Validated longitudinal-petiole phloem marker.",
  "Petiole longitudinal", "Phloem", "SEOR1", "PtXaAlbH.01G281100.v5.1", "sieve element", "main", "Validated longitudinal-petiole phloem marker.",
  "Petiole longitudinal", "Pith", "KT2", "PtXaTreH.15G030300.v5.1", "pith development", "main", "Validated longitudinal-petiole pith marker.",
  "Petiole longitudinal", "Pith", "STM", "PtXaAlbH.04G103400.v5.1", "meristematic identity", "main", "Validated longitudinal-petiole pith marker.",
  "Petiole longitudinal", "Vasculature1", "PXY", "PtXaAlbH.03G082400.v5.1", "procambium/cambium identity", "main", "Validated longitudinal-petiole Vasculature1 marker.",
  "Petiole longitudinal", "Vasculature1", "WOX4", "PtXaAlbH.14G016600.v5.1", "vascular cambium", "main", "Validated longitudinal-petiole Vasculature1 marker.",
  "Petiole longitudinal", "Vasculature2", "KAT1", "PtXaAlbH.04G003300.v5.1", "vascular transport", "main", "Validated longitudinal-petiole Vasculature2 marker.",
  "Petiole longitudinal", "Vasculature2", "UMAMIT19", "PtXaTreH.02G072800.v5.1", "amino-acid transport", "main", "Validated longitudinal-petiole Vasculature2 marker.",
  "Petiole longitudinal", "Xylem1", "CesA4", "PtXaTreH.02G213300.v5.1", "secondary cell wall", "main", "Validated longitudinal-petiole Xylem1 marker.",
  "Petiole longitudinal", "Xylem1", "LAC4", "PtXaAlbH.06G084500.v5.1", "xylem/lignification", "main", "Validated longitudinal-petiole Xylem1 marker.",
  "Petiole longitudinal", "Xylem2", "MYB59", "PtXaAlbH.01G196500.v5.1", "xylem differentiation", "main", "Validated longitudinal-petiole Xylem2 marker.",

  "Stem", "Cortex", "PNP", "PtXaAlbH.13G084300.v5.1", "domain-enriched", "main", "Dataset-supported cortex marker.",
  "Stem", "Cortex", "PMEAMT", "PtXaAlbH.15G031100.v5.1", "domain-enriched", "main", "Dataset-supported cortex marker.",
  "Stem", "Epidermis", "FDH", "PtXaTreH.06G179500.v5.1", "epidermis/cuticle", "main", "Strong epidermal/cuticle marker.",
  "Stem", "Epidermis", "KCS2", "PtXaAlbH.10G060600.v5.1", "epidermis/cuticle", "main", "Good epidermal/cuticle marker.",
  "Stem", "Pith", "WRKY12", "PtXaTreH.14G037200.v5.1", "pith/domain support", "main", "Preferred over HYDROLASE/PSBXb if detected.",
  "Stem", "Pith", "ATHB13", "PtXaTreH.10G076500.v5.1", "pith/domain support", "main", "Use if spatially enriched in pith.",
  "Stem", "Cambium", "CLE41", "PtXaTreH.02G198800.v5.1", "cambium/procambium identity", "main", "Preferred cambium-supporting marker.",
  "Stem", "Cambium", "ACL5a", "PtXaAlbH.06G184800.v5.1", "cambium/development", "main", "Use with CLE41 for cambium support.",
  "Stem", "Phloem", "SUT4", "PtXaTreH.02G091700.v5.1", "phloem/sucrose transport", "main", "Preferred phloem-supporting marker if detected.",
  "Stem", "Phloem", "UMAMIT12", "PtXaAlbH.06G071600.v5.1", "phloem/companion support", "main", "Added as second phloem-supporting marker if present.",
  "Stem", "Xylem", "XCP1", "PtXaAlbH.04G160300.v5.1", "xylem/vessel differentiation", "main", "Strong xylem/vessel marker.",
  "Stem", "Xylem", "LAC4", "PtXaAlbH.06G084500.v5.1", "xylem/lignification", "main", "Strong xylem/lignification marker.",
  "Stem", "Xylem", "4CL", "PtXaAlbH.01G031600.v5.1", "xylem/lignification", "supplement", "Good secondary-wall support marker."
)

dir.create("marker", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(fig2c_reviewed_markers, "marker/Fig2C_reviewed_curated_markers.csv")

audit_reviewed_markers <- function(marker_df, source_df) {
  if (is.null(source_df)) {
    return(marker_df %>% mutate(in_updated_markerlist = NA))
  }

  source_chr <- source_df %>%
    mutate(across(everything(), as.character)) %>%
    unlist(use.names = FALSE)

  marker_df %>%
    rowwise() %>%
    mutate(
      in_updated_markerlist = gene_id %in% source_chr ||
        marker_name %in% source_chr
    ) %>%
    ungroup()
}

fig2c_reviewed_markers_audited <- audit_reviewed_markers(fig2c_reviewed_markers, updated_marker_source)
readr::write_csv(fig2c_reviewed_markers_audited, "marker/Fig2C_reviewed_curated_markers_audited.csv")

plot_curated_marker_dotplot <- function(sobj, marker_df, out_pdf, group_col = "celltypes_v2",
                                        width = 7, height = 7, dot.scale = 6) {
  if (is.null(sobj)) {
    message("Reviewed marker dotplot skipped for ", out_pdf, ": Seurat object is NULL.")
    return(NULL)
  }
  stopifnot(group_col %in% colnames(sobj@meta.data))
  
  marker_df <- marker_df %>%
    filter(plot_set == "main") %>%
    filter(gene_id %in% rownames(sobj)) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  if (nrow(marker_df) == 0) stop("No marker genes found in the Seurat object.")
  
  Idents(sobj) <- group_col
  
  gene_labels <- setNames(paste0(marker_df$marker_name, "\n", marker_df$gene_id), marker_df$gene_id)
  
  p <- DotPlot_scCustom(
    seurat_object = sobj,
    features = marker_df$gene_id,
    group.by = group_col,
    flip_axes = TRUE,
    scale.by = "size",
    dot.min = 0,
    dot.scale = dot.scale,
    x_lab_rotate = FALSE
  ) +
    scale_y_discrete(labels = gene_labels) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 5.5, lineheight = 0.9),
      axis.title = element_blank(),
      legend.position = "right"
    ) +
    labs(color = "Average expression", size = "Percent expressed")
  
  ggsave(out_pdf, p, width = width, height = height)
}

plot_curated_marker_heatmap <- function(sobj, marker_df, out_pdf, group_col = "celltypes_v2",
                                        width = 7, height = 7, assay_use = NULL) {
  if (is.null(sobj)) {
    message("Reviewed marker heatmap skipped for ", out_pdf, ": Seurat object is NULL.")
    return(NULL)
  }
  stopifnot(group_col %in% colnames(sobj@meta.data))

  marker_df <- marker_df %>%
    filter(plot_set == "main") %>%
    filter(gene_id %in% rownames(sobj)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    mutate(gene_label = paste0(marker_name, "\n", gene_id))

  if (nrow(marker_df) == 0) stop("No marker genes found in the Seurat object.")

  if (is.null(assay_use)) {
    assay_use <- if ("SCT" %in% names(sobj@assays)) "SCT" else DefaultAssay(sobj)
  }
  avg_mtx <- AverageExpression(sobj, assay = assay_use, group.by = group_col, verbose = FALSE)[[assay_use]]
  avg_mtx <- avg_mtx[marker_df$gene_id, , drop = FALSE]
  z_mtx <- t(scale(t(avg_mtx)))
  z_mtx[!is.finite(z_mtx)] <- 0
  z_mtx <- pmax(pmin(z_mtx, 2), -2)

  # [2026.7.24 FIX] check.names = FALSE: domain names contain spaces
  #   ("Shoot meristem"); base as.data.frame() mangles them to "Shoot.meristem",
  #   which then fails to match the factor levels -> all-NA (blank) heatmap.
  plot_df <- as.data.frame(z_mtx, check.names = FALSE) %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "celltype", values_to = "z_score") %>%
    left_join(marker_df %>% select(gene_id, gene_label), by = "gene_id") %>%
    mutate(
      gene_label = factor(gene_label, levels = rev(marker_df$gene_label)),
      celltype = factor(celltype, levels = colnames(z_mtx))
    )

  p <- ggplot(plot_df, aes(x = celltype, y = gene_label, fill = z_score)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-2, 2), name = "Scaled\naverage\nexpression"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 5.5, lineheight = 0.9),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    )

  ggsave(out_pdf, p, width = width, height = height)
}

sam_curated_markers <- fig2c_reviewed_markers_audited %>% filter(tissue == "SAM")
bud_curated_markers <- fig2c_reviewed_markers_audited %>% filter(tissue == "Axillary bud")
petiole_curated_markers <- fig2c_reviewed_markers_audited %>% filter(tissue == "Petiole cross")
petiole_long_curated_markers <- fig2c_reviewed_markers_audited %>% filter(tissue == "Petiole longitudinal")
stem_curated_markers <- fig2c_reviewed_markers_audited %>% filter(tissue == "Stem")

plot_curated_marker_dotplot(sobj_sam_ann, sam_curated_markers, "summary/Fig2_annotation/Fig2C_SAM_reviewed_marker_dotplot.pdf", group_col = "celltype", width = 7, height = 8)
plot_curated_marker_dotplot(sobj_bud_ann, bud_curated_markers, "summary/Fig2_annotation/Fig2C_bud_reviewed_marker_dotplot.pdf", group_col = "celltype", width = 7, height = 8)
plot_curated_marker_dotplot(sobj_petiole_cross_ann, petiole_curated_markers, "summary/Fig2_annotation/Fig2C_petiole_cross_reviewed_marker_dotplot.pdf", group_col = "celltype", width = 6.5, height = 5.5)
plot_curated_marker_dotplot(sobj_stem_cross_ann, stem_curated_markers, "summary/Fig2_annotation/Fig2C_stem_reviewed_marker_dotplot.pdf", group_col = "celltype", width = 7, height = 7)

plot_curated_marker_heatmap(sobj_sam_ann, sam_curated_markers, "summary/Fig2_annotation/Fig2C_SAM_reviewed_marker_heatmap.pdf", group_col = "celltype", width = 7, height = 8)
plot_curated_marker_heatmap(sobj_bud_ann, bud_curated_markers, "summary/Fig2_annotation/Fig2C_bud_reviewed_marker_heatmap.pdf", group_col = "celltype", width = 7, height = 8)
plot_curated_marker_heatmap(sobj_petiole_cross_ann, petiole_curated_markers, "summary/Fig2_annotation/Fig2C_petiole_cross_reviewed_marker_heatmap.pdf", group_col = "celltype", width = 6.5, height = 5.5)
plot_curated_marker_heatmap(sobj_petiole_long_ann, petiole_long_curated_markers, "summary/Fig2_annotation/Fig2C_petiole_long_reviewed_marker_heatmap.pdf", group_col = "celltype", width = 7, height = 6.5)
plot_curated_marker_heatmap(sobj_stem_cross_ann, stem_curated_markers, "summary/Fig2_annotation/Fig2C_stem_reviewed_marker_heatmap.pdf", group_col = "celltype", width = 7, height = 7)
# MANUAL STITCH for Fig. 2:

#   Combine Fig2A spatial maps + Fig2B UMAPs + representative Fig2C dotplots.
#   Keep full dotplots for all tissues in supplementary if too large.

# --------------------------------------- #
# 3. Marker atlas figure ####
# --------------------------------------- #
#
# Fig. 3 tells one story in four steps:
#   3A  known markers from the literature are domain-specific across major tissues
#   3B  three of those known markers resolve in situ in stem cross sections
#   3C  de novo markers for the anchor tissue (stem), one dotplot
#   3D  four de novo markers mapped in situ on the same stem sections
#
# Design rules that constrain the code below:
#   - a marker "works" if it is specific to its expected domain within a tissue;
#     off-target expression in an unrelated domain does not disqualify it
#   - tissues are annotated at different resolutions, so marker x tissue cells
#     that cannot be tested are drawn as grey states, never as a value
#   - remaining tissues and the full 29-marker audit go to supplementary
#
# Everything below runs from the workspace root. Panels go to
# summary/Fig3_marker_atlas/, supplements to summary/Supplementary/FigS_marker_atlas/,
# and expensive intermediates are cached in temporary/fig3_cache/.

FIG3_DIR  <- "summary/Fig3_marker_atlas"
FIG3_SUPP <- "summary/Supplementary/FigS_marker_atlas"
FIG3_CACHE <- "temporary/fig3_cache"
for (d in c(FIG3_DIR, FIG3_SUPP, FIG3_CACHE)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---------- 3.0 shared vocabulary and helpers ----------

# The five tissues were annotated independently and use different label strings
# for the same biology. Fig. 3 needs one axis, so every label is mapped onto a
# single controlled vocabulary. Anything unmapped becomes NA and is dropped:
# silently lumping unknown labels into a domain would corrupt every panel.
FIG3_DOMAIN_LEVELS <- c(
  "Epidermis", "Cortex", "Inner cortex", "Pith", "Pith (Xylem)", "Mesophyll",
  "Bud scale", "Leaf primordium", "Bud primordium", "Proliferating",
  "Meristem", "Shoot meristem", "Cambium", "Bud procambium", "Vasculature",
  "Vasculature (Phloem)", "Phloem", "Xylem"
)

harmonize_domain <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("Epidermis", "Epidermis1", "Epidermis2", "Dividing epidermis")        ~ "Epidermis",
    x %in% c("Cortex", "Cortex1", "Cortex2", "Inner cortex", "Outer cortex",
             "Cortex/ground parenchyma", "Ground tissue",
             "Collenchyma", "Parenchyma")                                          ~ "Cortex",
    x %in% c("Pith", "Pith parenchyma")                                            ~ "Pith",
    x %in% c("Mesophyll", "Mature mesophyll")                                      ~ "Mesophyll",
    x %in% c("Bud scale", "Bud scales")                                            ~ "Bud scale",
    x %in% c("Leaf primordium", "Bud primordium", "Primordium",
             "Leaf/bud primordium")                                                ~ "Primordium",
    x %in% c("Proliferating", "Proliferating cells")                               ~ "Proliferating",
    x %in% c("Meristem", "Meristematic", "SAM", "Shoot meristem",
             "Shoot apical meristem")                                              ~ "Meristem",
    x %in% c("Cambium", "Procambium", "Cambium/procambium", "Bud procambium",
             "Vascular cambium", "Procambial stem cells")                          ~ "Cambium/procambium",
    x %in% c("Vasculature", "Vasculature1", "Vasculature2",
             "Vascular bundle", "Vasculature (mixed)")                             ~ "Vasculature (mixed)",
    x %in% c("Phloem", "Vasculature (Phloem)")                                     ~ "Phloem",
    x %in% c("Xylem", "Xylem1", "Xylem2", "Pith (Xylem)")                          ~ "Xylem",
    TRUE                                                                           ~ NA_character_
  )
}

strip_ver <- function(x) sub("\\.v5\\.1$", "", x)

FIG3_BASE <- 7; FIG3_SMALL <- 6; FIG3_TICK <- 5.5
theme_fig3 <- function() {
  theme_bw(base_size = FIG3_BASE) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = FIG3_TICK, colour = "black"),
      axis.title = element_text(size = FIG3_BASE),
      strip.background = element_rect(fill = "grey94", colour = NA),
      strip.text = element_text(size = FIG3_SMALL),
      legend.key.size = unit(3, "mm"),
      plot.title = element_text(size = FIG3_BASE, face = "plain")
    )
}

FIG3_EXPR_COLS <- c("#4575B4", "#E8EDF3", "#FDDBC7", "#F4A582", "#B2182B")
scale_fill_fig3_z <- function(lo = -1.4, hi = 2.2) {
  scale_fill_gradientn(
    colours = FIG3_EXPR_COLS,
    values = scales::rescale(c(lo, lo * 0.25, 0.35, hi * 0.55, hi)),
    limits = c(lo, hi), breaks = c(-1, 0, 1, 2), name = "Expression (z)",
    guide = guide_colourbar(barheight = unit(2.2, "mm"), barwidth = unit(14, "mm"),
                            title.position = "top", ticks = FALSE)
  )
}

# Per-domain mean expression, detection fraction, and a domain-wise z-score.
# Computed once per tissue and cached: the objects are large and every panel
# needs these three matrices.
FIG3_CACHE_VERSION <- "final_celltype_labels_v1"

build_fig3_cache <- function(obj, tissue, key, anno_col) {
  f <- file.path(FIG3_CACHE, paste0("cache_", key, ".rds"))
  if (file.exists(f)) {
    old <- readRDS(f)
    if (isTRUE(old$cache_version == FIG3_CACHE_VERSION) && identical(old$anno_col, anno_col)) return(old)
  }
  DefaultAssay(obj) <- "SCT"
  dom <- as.character(obj@meta.data[[anno_col]])
  keep <- !is.na(dom) & nzchar(dom) & dom %in% FIG3_DOMAIN_LEVELS
  m <- GetAssayData(obj, assay = "SCT", layer = "data")[, keep, drop = FALSE]
  dom <- dom[keep]
  doms <- intersect(FIG3_DOMAIN_LEVELS, unique(dom))
  avg <- sapply(doms, function(d) Matrix::rowMeans(m[, dom == d, drop = FALSE]))
  pct <- sapply(doms, function(d) Matrix::rowMeans(m[, dom == d, drop = FALSE] > 0))
  out <- list(tissue = tissue, key = key, anno_col = anno_col,
              cache_version = FIG3_CACHE_VERSION,
              avg = avg, pct = pct,
              n_spots = table(factor(dom, levels = doms)),
              genes = rownames(m), images = Images(obj))
  saveRDS(out, f)
  out
}

fig3_z <- function(cache) {
  z <- t(scale(t(as.matrix(cache$avg))))
  z[!is.finite(z)] <- 0
  z
}

# Specificity: how much higher a gene is in its intended domain than in the
# best competing domain. This is the filter that stops "expressed everywhere,
# significant by DE" genes from being called markers.
fig3_specificity <- function(gene, cache, intended) {
  a <- cache$avg; p <- cache$pct
  doms <- intersect(FIG3_DOMAIN_LEVELS, colnames(a))
  ii <- intersect(intended, doms)
  if (!gene %in% rownames(a) || !length(ii)) {
    return(tibble(tgt_expr = NA_real_, tgt_pct = NA_real_,
                  max_off_expr = NA_real_, max_off_pct = NA_real_, spec_log2fc = NA_real_))
  }
  off <- setdiff(doms, ii)
  te <- max(a[gene, ii]); oe <- if (length(off)) max(a[gene, off]) else 0
  tibble(tgt_expr = te, tgt_pct = max(p[gene, ii]),
         max_off_expr = oe, max_off_pct = if (length(off)) max(p[gene, off]) else 0,
         spec_log2fc = log2((te + 0.05) / (oe + 0.05)))
}

# Display identifiers: strip the haplotype prefix and the annotation version.
# Both haplotypes appear in the marker sets (PtXaAlbH., PtXaTreH.), so the
# pattern has to match either -- anchoring on one leaves half the axis labels
# carrying a prefix the rest do not have.
fig3_short_id <- function(x) sub("\\.v5\\.1$", "", sub("^Pt[A-Za-z]+\\.", "", x))

FIG3_OBJ_SPEC <- list(
  sam           = list(tissue = "SAM",                  anno = "celltype"),
  bud           = list(tissue = "Axillary bud",         anno = "celltype"),
  stem_cross    = list(tissue = "Stem",                 anno = "celltype"),
  petiole_cross = list(tissue = "Petiole cross",        anno = "celltype"),
  petiole_long  = list(tissue = "Petiole longitudinal", anno = "celltype")
)

fig3_objs <- list(
  sam           = sobj_sam_ann,
  bud           = sobj_bud_ann,
  stem_cross    = sobj_stem_cross_ann,
  petiole_cross = sobj_petiole_cross_ann,
  petiole_long  = sobj_petiole_long_ann
)
for (k in names(fig3_objs)) {
  fig3_objs[[k]]$fig3_domain <- as.character(fig3_objs[[k]]@meta.data[[FIG3_OBJ_SPEC[[k]]$anno]])
  fig3_domain_levels <- intersect(FIG3_DOMAIN_LEVELS, unique(as.character(fig3_objs[[k]]$fig3_domain)))
  fig3_objs[[k]]$fig3_domain <- factor(fig3_objs[[k]]$fig3_domain, levels = fig3_domain_levels)
}

fig3_caches <- lapply(names(FIG3_OBJ_SPEC), function(k)
  build_fig3_cache(fig3_objs[[k]], FIG3_OBJ_SPEC[[k]]$tissue, k, FIG3_OBJ_SPEC[[k]]$anno))
names(fig3_caches) <- names(FIG3_OBJ_SPEC)
fig3_z_list <- lapply(fig3_caches, fig3_z)

# Every original annotation label, what it maps to, and how many spots it
# carries. This is the audit trail for the harmonization: a reader can check
# that e.g. "Pith (Xylem)" in petiole cross was mapped to Xylem deliberately.
fig3_harmonization_map <- purrr::map_dfr(names(FIG3_OBJ_SPEC), function(k) {
  ac <- FIG3_OBJ_SPEC[[k]]$anno
  lab <- as.character(fig3_objs[[k]]@meta.data[[ac]])
  tibble(tissue = FIG3_OBJ_SPEC[[k]]$tissue, annotation_column = ac,
         original_label = lab) %>%
    count(tissue, annotation_column, original_label, name = "n_spots") %>%
    mutate(harmonized_domain = harmonize_domain(original_label)) %>%
    arrange(harmonized_domain, original_label) %>%
    select(tissue, annotation_column, original_label, harmonized_domain, n_spots)
})
readr::write_csv(fig3_harmonization_map,
                 file.path(FIG3_DIR, "Fig3_domain_harmonization_map.csv"))
stopifnot(!any(is.na(fig3_harmonization_map$harmonized_domain)))   # unmapped label = silent data loss


# The curated known-marker set. known_marker_audit.rds is built by the
# marker-curation block: it assembles the four project sources, resolves
# synonyms, assigns Arabidopsis orthologs, and attaches each marker's intended
# domain(s) in `domain_accept` (multiple domains are PIPE-separated -- splitting
# on the wrong character silently reduces every multi-domain marker to one
# unmatched target). `main_sel` is the 29-marker set that survived curation.
fig3_marker_audit <- readRDS(file.path(FIG3_CACHE, "known_marker_audit.rds"))
main_sel <- fig3_marker_audit$main_sel
stopifnot(all(c("gene_id", "marker_name", "domain_accept", "marker_class",
                "evidence_source", "source_detail", "ath_hit", "func_anno")
              %in% names(main_sel)))

FIG3_DE_TABLES <- c(
  sam = NA_character_,
  bud = "marker/BUD_all_markers_by_celltype_v2_2026_6.csv",
  stem_cross = "marker/STEM_cross_AB_all_markers_by_celltype_PtXaOnly.csv",
  petiole_cross = "marker/PETIOLE_CROSS_all_markers_by_celltype_PtXaOnly.csv",
  petiole_long = "marker/PETIOLE_LONG_all_markers_by_celltype_PtXaOnly.csv"
)

# Per-domain differential expression, cached per tissue. Domains come from the
# final celltype labels so Fig. 3 matches the Fig. 2 annotation vocabulary.
fig3_run_de <- function(obj, key, anno_col, prep_sct = FALSE) {
  f <- file.path(FIG3_CACHE, paste0("de_", key, ".rds"))
  external_file <- FIG3_DE_TABLES[[key]]
  if (!is.na(external_file) && file.exists(external_file)) {
    d <- read.csv(external_file, check.names = FALSE)
    if (!"gene_id" %in% names(d) && "gene" %in% names(d)) d$gene_id <- d$gene
    d <- d %>% filter(cluster %in% unique(as.character(obj@meta.data[[anno_col]])))
    attr(d, "cache_version") <- FIG3_CACHE_VERSION
    saveRDS(d, f)
    return(d)
  }
  if (file.exists(f)) {
    old <- readRDS(f)
    if (isTRUE(attr(old, "cache_version") == FIG3_CACHE_VERSION) &&
        "cluster" %in% names(old) &&
        all(c("gene_id", "p_val_adj", "avg_log2FC", "pct.1", "pct.2") %in% names(old)) &&
        all(unique(as.character(old$cluster)) %in% FIG3_DOMAIN_LEVELS)) return(old)
  }
  obj$fig3_domain <- as.character(obj@meta.data[[anno_col]])
  obj <- obj[, !is.na(obj$fig3_domain) & obj$fig3_domain %in% FIG3_DOMAIN_LEVELS]
  DefaultAssay(obj) <- "SCT"
  future::plan("sequential")
  options(future.globals.maxSize = 32 * 1024^3)
  if (prep_sct) {
    obj <- PrepSCTFindMarkers(obj, verbose = FALSE)
  }
  Idents(obj) <- "fig3_domain"
  d <- FindAllMarkers(obj, assay = "SCT", slot = "data", only.pos = TRUE,
                      min.pct = 0.10, logfc.threshold = 0.25,
                      recorrect_umi = FALSE, verbose = FALSE)
  if (nrow(d) == 0) {
    d <- tibble(gene_id = character(), cluster = character(), p_val_adj = numeric(),
                avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric())
    attr(d, "cache_version") <- FIG3_CACHE_VERSION
    saveRDS(d, f)
    return(d)
  }
  if (!"gene" %in% names(d)) d <- tibble::rownames_to_column(d, "gene")
  d <- d %>% rename(gene_id = gene)
  attr(d, "cache_version") <- FIG3_CACHE_VERSION
  saveRDS(d, f)
  d
}

# ---------- 3.1 Fig. 3A: published markers across major tissues ----------
#
# Fig. 3 uses a compact main story. Published markers are shown across the major
# annotated tissues to demonstrate that conserved markers peak in the expected
# domains. For coarse annotations, vascular markers are allowed to peak in a
# broader vascular compartment; the full marker audit is supplementary.
#
#   margin = log2( (max mean expr over accepted target domains + 0.05) /
#                  (max mean expr over all other domains       + 0.05) )
#
# Three states are distinguished and drawn differently, because collapsing them
# would overstate coverage:
#   tested               - the intended domain exists here and the gene is present
#   domain_not_annotated - this tissue has no such domain (light grey)
#   gene_absent          - the gene is not in this tissue's matrix (dark grey)
#
# "validated" (white dot) = margin > 0 AND >=10% of target-domain spots detect
# the gene. Off-target expression elsewhere is not penalised: the claim is
# within-tissue domain specificity, not plant-wide uniqueness.

fig3_expand_accept_domains <- function(accept_domains) {
  x <- trimws(accept_domains)
  expanded <- purrr::map(x, function(d) {
    switch(d,
           "Epidermis" = c("Epidermis"),
           "Cortex" = c("Cortex", "Inner cortex"),
           "Pith" = c("Pith", "Pith (Xylem)"),
           "Mesophyll" = c("Mesophyll"),
           "Bud scale" = c("Bud scale"),
           "Primordium" = c("Leaf primordium", "Bud primordium"),
           "Proliferating" = c("Proliferating"),
           "Meristem" = c("Meristem", "Shoot meristem"),
           "Cambium/procambium" = c("Cambium", "Bud procambium"),
           "Vasculature (mixed)" = c("Vasculature", "Vasculature (Phloem)"),
           "Phloem" = c("Phloem", "Vasculature (Phloem)"),
           "Xylem" = c("Xylem", "Pith (Xylem)"),
           d)
  }) %>% unlist(use.names = FALSE)
  unique(c(x, expanded))
}

fig3_spec_margin <- function(gene_id, cache, accept_domains) {
  doms <- colnames(cache$avg)
  tg   <- intersect(fig3_expand_accept_domains(accept_domains), doms)
  if (!length(tg))
    return(tibble(status = "domain_not_annotated", margin = NA_real_,
                  best_target = NA_character_, peak_domain = NA_character_,
                  pct_target = NA_real_))
  if (!gene_id %in% rownames(cache$avg))
    return(tibble(status = "gene_absent", margin = NA_real_,
                  best_target = NA_character_, peak_domain = NA_character_,
                  pct_target = NA_real_))
  v   <- cache$avg[gene_id, ]
  p   <- cache$pct[gene_id, ]
  off <- setdiff(doms, tg)
  bt  <- tg[which.max(v[tg])]
  tibble(status      = "tested",
         margin      = log2((max(v[tg]) + 0.05) /
                            (if (length(off)) max(v[off]) + 0.05 else 0.05)),
         best_target = bt,
         peak_domain = doms[which.max(v)],
         pct_target  = p[bt])
}

# audit every curated marker in every tissue (this is also the supplementary table)
fig3_known_audit <- purrr::map_dfr(names(FIG3_OBJ_SPEC), function(k) {
  purrr::pmap_dfr(
    list(main_sel$gene_id, main_sel$marker_name,
         main_sel$domain_accept, main_sel$marker_class),
    function(gi, mn, acc, cls) {
      bind_cols(
        tibble(marker_name = mn, gene_id = gi, marker_class = cls,
               intended_domain = acc, tissue = FIG3_OBJ_SPEC[[k]]$tissue),
        fig3_spec_margin(gi, fig3_caches[[k]], trimws(strsplit(acc, "\\|")[[1]])))
    })
}) %>%
  mutate(validated = status == "tested" & margin > 0 & pct_target >= 0.10)

# compress multi-domain targets to a short row label, e.g. "Cortex|Pith" -> "Cortex/pith"
fig3_short_target <- function(x) {
  v <- trimws(strsplit(x, "\\|")[[1]])
  if (length(v) == 1) v else paste0(v[1], "/", tolower(sub(" .*", "", v[2])))
}

FIG3A_CLASS_ORDER <- c("Epidermis", "Ground tissue", "Meristem / primordium",
                       "Cambium / procambium", "Phloem", "Xylem")
FIG3A_TISSUES <- c("SAM", "Axillary bud", "Stem", "Petiole cross")
FIG3A_TISSUE_COLS <- c("SAM" = "#E76BF3", "Axillary bud" = "#F8766D",
                       "Stem" = "#A3A500", "Petiole cross" = "#00B0F4")
FIG3A_MARKER_GROUPS <- tibble::tribble(
  ~marker_name, ~marker_type,
  "KCS2",       "Epidermis",
  "FDH",        "Epidermis",
  "LPTG1",      "Epidermis",
  "PDF1",       "Epidermis",
  "WRKY12",     "Ground tissue",
  "ATHB13",     "Ground tissue",
  "ANTb",       "Meristem / primordium",
  "CYCB1;5",    "Meristem / primordium",
  "MCM7",       "Meristem / primordium",
  "PXY",        "Cambium / procambium",
  "WOX4",       "Cambium / procambium",
  "CLE41",      "Cambium / procambium",
  "AtHB8",      "Cambium / procambium",
  "PP2-A10-2",  "Phloem",
  "SEOR1",      "Phloem",
  "APL",        "Phloem",
  "4CL",        "Xylem",
  "XCP2",       "Xylem",
  "MAN6a/b",    "Xylem",
  "LAC4",       "Xylem",
  "CESA8",      "Xylem"
)

fig3a_compatible_domains <- function(marker_type, domain_accept) {
  acc <- trimws(strsplit(domain_accept, "\\|")[[1]])
  extra <- dplyr::case_when(
    marker_type %in% c("Cambium / procambium", "Phloem", "Xylem") ~ "Vasculature (mixed)",
    TRUE ~ NA_character_
  )
  fig3_expand_accept_domains(unique(na.omit(c(acc, extra))))
}

fig3a_candidate_markers <- main_sel %>%
  inner_join(FIG3A_MARKER_GROUPS, by = "marker_name") %>%
  mutate(marker_type = factor(marker_type, levels = FIG3A_CLASS_ORDER)) %>%
  arrange(marker_type, marker_name)

fig3a_marker_rank <- purrr::pmap_dfr(
  list(fig3a_candidate_markers$gene_id, fig3a_candidate_markers$marker_name,
       fig3a_candidate_markers$domain_accept, fig3a_candidate_markers$marker_class,
       fig3a_candidate_markers$marker_type),
  function(gi, mn, acc, cls, mt) {
    acc2 <- fig3a_compatible_domains(as.character(mt), acc)
    tissue_stats <- purrr::map_dfr(names(FIG3_OBJ_SPEC), function(k) {
      tissue <- FIG3_OBJ_SPEC[[k]]$tissue
      if (!tissue %in% FIG3A_TISSUES) return(NULL)
      fig3_spec_margin(gi, fig3_caches[[k]], acc2) %>% mutate(tissue = tissue)
    })
    tibble(marker_name = mn, gene_id = gi, marker_type = as.character(mt), marker_class = cls,
           domain_accept = acc, compatible_domain = paste(acc2, collapse = "|"),
           n_tested = sum(tissue_stats$status == "tested"),
           n_supported = sum(tissue_stats$status == "tested" & tissue_stats$margin > 0 & tissue_stats$pct_target >= 0.10),
           median_margin = median(tissue_stats$margin[tissue_stats$status == "tested"], na.rm = TRUE),
           mean_pct_target = mean(tissue_stats$pct_target[tissue_stats$status == "tested"], na.rm = TRUE),
           evidence_source = main_sel$evidence_source[main_sel$gene_id == gi][1],
           source_detail = main_sel$source_detail[main_sel$gene_id == gi][1],
           arabidopsis_best_hit = main_sel$ath_hit[main_sel$gene_id == gi][1],
           functional_annotation = main_sel$func_anno[main_sel$gene_id == gi][1])
  }) %>%
  mutate(median_margin = ifelse(is.finite(median_margin), median_margin, NA_real_),
         mean_pct_target = ifelse(is.finite(mean_pct_target), mean_pct_target, NA_real_),
         marker_type = factor(marker_type, levels = FIG3A_CLASS_ORDER),
         conservation_score = n_supported * 100 + n_tested * 10 +
           coalesce(median_margin, 0) + coalesce(mean_pct_target, 0)) %>%
  arrange(marker_type, desc(conservation_score), desc(n_supported), desc(median_margin)) %>%
  group_by(marker_type) %>%
  mutate(class_rank = row_number()) %>%
  ungroup()

FIG3A_MARKERS_TO_PLOT <- c(
  "KCS2", "FDH",
  "WRKY12", "ATHB13",
  "ANTb", "CYCB1;5",
  "PXY", "AtHB8",
  "PP2-A10-2", "SEOR1",
  "4CL", "MAN6a/b"
)
fig3a_known_sel <- fig3a_marker_rank %>%
  filter(marker_name %in% FIG3A_MARKERS_TO_PLOT, n_supported >= 2) %>%
  mutate(plot_rank = match(marker_name, FIG3A_MARKERS_TO_PLOT)) %>%
  arrange(plot_rank)
fig3a_keep <- as.character(fig3a_known_sel$marker_name)
readr::write_csv(
  fig3a_marker_rank %>%
    mutate(marker_type = as.character(marker_type),
           marker_class = as.character(marker_class),
           plotted_in_Fig3A = marker_name %in% fig3a_keep) %>%
    select(marker_name, gene_id, marker_type, marker_class, domain_accept, compatible_domain,
           n_tested, n_supported, median_margin, mean_pct_target,
           conservation_score, class_rank, plotted_in_Fig3A,
           evidence_source, source_detail, arabidopsis_best_hit, functional_annotation),
  file.path(FIG3_DIR, "Fig3A_published_marker_conservation_ranked_table.csv"))

readr::write_csv(
  fig3a_marker_rank %>%
    filter(marker_type %in% c("Cambium / procambium", "Phloem", "Xylem")) %>%
    mutate(marker_type = as.character(marker_type),
           marker_class = as.character(marker_class),
           plotted_in_Fig3A = marker_name %in% fig3a_keep) %>%
    select(marker_name, gene_id, marker_type, marker_class, compatible_domain,
           n_tested, n_supported, median_margin, mean_pct_target,
           conservation_score, class_rank, plotted_in_Fig3A,
           evidence_source, source_detail, arabidopsis_best_hit, functional_annotation),
  file.path(FIG3_DIR, "Fig3A_known_vascular_xylem_candidate_audit.csv"))

fig3a_plot_df <- purrr::map_dfr(names(FIG3_OBJ_SPEC), function(k) {
  cache <- fig3_caches[[k]]
  tissue <- FIG3_OBJ_SPEC[[k]]$tissue
  if (!tissue %in% FIG3A_TISSUES) return(NULL)
  purrr::map_dfr(seq_len(nrow(fig3a_known_sel)), function(i) {
    gi <- fig3a_known_sel$gene_id[i]
    mn <- as.character(fig3a_known_sel$marker_name[i])
    mt <- as.character(fig3a_known_sel$marker_type[i])
    cls <- as.character(fig3a_known_sel$marker_class[i])
    acc <- trimws(strsplit(fig3a_known_sel$compatible_domain[i], "\\|")[[1]])
    tibble(
      tissue = tissue,
      domain = colnames(cache$avg),
      marker_name = mn,
      gene_id = gi,
      marker_type = mt,
      marker_class = cls,
      expected = domain %in% acc,
      avg_expr = if (gi %in% rownames(cache$avg)) as.numeric(cache$avg[gi, ]) else 0,
      pct_exp = if (gi %in% rownames(cache$pct)) 100 * as.numeric(cache$pct[gi, ]) else 0
    )
  })
}) %>%
  mutate(tissue = factor(tissue, levels = FIG3A_TISSUES),
         domain = factor(domain, levels = FIG3_DOMAIN_LEVELS),
         x_lab = paste(tissue, domain, sep = " | "),
         y_lab = paste0(marker_name, "\n", gene_id))

fig3a_y_levels <- fig3a_known_sel %>%
  mutate(y_lab = paste0(marker_name, "\n", gene_id)) %>%
  pull(y_lab) %>%
  rev()
make_fig3a_marker_dotplot <- function(plot_df, tissues, file_prefix,
                                      width = 6, height = 5,
                                      flip_width = 7, flip_height = 4.5) {
  d <- plot_df %>%
    filter(as.character(tissue) %in% tissues) %>%
    mutate(tissue = factor(as.character(tissue), levels = tissues))
  x_levels <- d %>%
    distinct(tissue, domain, x_lab) %>%
    arrange(tissue, domain) %>%
    pull(x_lab)
  x_map <- tibble(x_lab = x_levels, x_pos = seq_along(x_levels),
                  x_short = sub("^.* \\| ", "", x_levels))
  d <- d %>%
    left_join(x_map, by = "x_lab") %>%
    group_by(gene_id) %>%
    mutate(avg_scaled = as.numeric(scale(avg_expr))) %>%
    ungroup() %>%
    mutate(avg_scaled = ifelse(is.finite(avg_scaled), avg_scaled, 0),
           avg_scaled = pmax(pmin(avg_scaled, 2), -2),
           y_lab = factor(y_lab, levels = fig3a_y_levels),
           y_pos = as.numeric(y_lab))
  tissue_bar <- d %>%
    distinct(tissue, x_pos) %>%
    group_by(tissue) %>%
    summarise(xmin = min(x_pos) - 0.5,
              xmax = max(x_pos) + 0.5,
              xmid = mean(range(x_pos)),
              .groups = "drop")
  tissue_boundaries <- tissue_bar %>%
    arrange(factor(tissue, levels = tissues))
  tissue_boundaries <- tissue_boundaries[
    seq_len(max(nrow(tissue_boundaries) - 1, 0)), , drop = FALSE]
  marker_bar <- d %>%
    distinct(y_lab, y_pos, marker_type) %>%
    group_by(marker_type) %>%
    summarise(ymin = min(y_pos) - 0.5,
              ymax = max(y_pos) + 0.5,
              ymid = mean(range(y_pos)),
              .groups = "drop") %>%
    mutate(marker_type_lab = gsub(" / ", "\n", as.character(marker_type)))

  p_tissue_bar <- ggplot(tissue_bar) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0.5, ymax = 1.5, fill = tissue),
              colour = "white", linewidth = 0.25) +
    geom_segment(data = tissue_boundaries,
                 aes(x = xmax, xend = xmax, y = 0.5, yend = 1.5),
                 inherit.aes = FALSE, colour = "grey35", linewidth = 0.3) +
    geom_text(aes(x = xmid, y = 1, label = tissue), size = 2.1, colour = "black") +
    scale_fill_manual(values = FIG3A_TISSUE_COLS, guide = "none") +
    scale_x_continuous(limits = c(0.5, length(x_levels) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = FIG3_SMALL) +
    theme(plot.margin = margin(0, 2, 0, 2))
  p_marker_bar <- ggplot(marker_bar) +
    geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = ymin, ymax = ymax),
              fill = "grey88", colour = "white", linewidth = 0.25) +
    geom_text(aes(x = 1, y = ymid, label = marker_type_lab),
              size = 1.8, colour = "black", lineheight = 0.85) +
    scale_y_continuous(limits = c(0.5, length(fig3a_y_levels) + 0.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = FIG3_SMALL) +
    theme(plot.margin = margin(0, 0, 2, 2))
  p_dot <- ggplot(d, aes(x_pos, y_pos)) +
    geom_vline(data = tissue_boundaries, aes(xintercept = xmax),
               inherit.aes = FALSE, colour = "grey65", linewidth = 0.25) +
    geom_point(aes(size = pct_exp, colour = avg_scaled), alpha = 0.95) +
    scale_colour_gradient2(low = "#3B4CC0", mid = "grey92", high = "#B40426",
                           midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                           name = "Scaled\nexpr.") +
    scale_size_area(max_size = 2.8, limits = c(0, 100), name = "% spots") +
    scale_x_continuous(breaks = x_map$x_pos, labels = x_map$x_short,
                       limits = c(0.5, length(x_levels) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(fig3a_y_levels), labels = fig3a_y_levels,
                       limits = c(0.5, length(fig3a_y_levels) + 0.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = FIG3_SMALL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.8),
          axis.text.y = element_text(size = 4.2, lineheight = 0.88),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5.2),
          legend.key.size = unit(7, "pt"),
          plot.margin = margin(0, 2, 2, 2))
  p_main <- (((patchwork::plot_spacer() | p_tissue_bar) +
                patchwork::plot_layout(widths = c(0.18, 1))) /
             ((p_marker_bar | p_dot) +
                patchwork::plot_layout(widths = c(0.18, 1)))) +
    patchwork::plot_layout(heights = c(0.08, 1), guides = "keep")
  ggsave(file.path(FIG3_DIR, paste0(file_prefix, ".pdf")), p_main,
         width = width, height = height, device = cairo_pdf)

  d_flip <- d %>%
    mutate(marker_lab = factor(as.character(y_lab), levels = rev(fig3a_y_levels)),
           marker_pos = as.numeric(marker_lab),
           domain_lab = factor(x_lab, levels = rev(x_levels)))
  flip_bar <- d_flip %>%
    distinct(marker_pos, marker_type) %>%
    group_by(marker_type) %>%
    summarise(xmin = min(marker_pos) - 0.5,
              xmax = max(marker_pos) + 0.5,
              xmid = mean(range(marker_pos)),
              .groups = "drop") %>%
    mutate(marker_type_lab = gsub(" / ", "\n", as.character(marker_type)))
  p_flip_bar <- ggplot(flip_bar) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0.5, ymax = 1.5),
              fill = "grey88", colour = "white", linewidth = 0.25) +
    geom_text(aes(x = xmid, y = 1, label = marker_type_lab),
              size = 1.8, colour = "black", lineheight = 0.85) +
    scale_x_continuous(limits = c(0.5, length(fig3a_y_levels) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = FIG3_SMALL) +
    theme(plot.margin = margin(0, 2, 0, 2))
  p_flip_dot <- ggplot(d_flip, aes(marker_pos, domain_lab)) +
    geom_point(aes(size = pct_exp, colour = avg_scaled), alpha = 0.95) +
    scale_colour_gradient2(low = "#3B4CC0", mid = "grey92", high = "#B40426",
                           midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                           name = "Scaled\nexpr.") +
    scale_size_area(max_size = 2.8, limits = c(0, 100), name = "% spots") +
    scale_x_continuous(breaks = seq_along(rev(fig3a_y_levels)),
                       labels = rev(fig3a_y_levels),
                       limits = c(0.5, length(fig3a_y_levels) + 0.5), expand = c(0, 0)) +
    scale_y_discrete(labels = setNames(paste0(sub(" \\| .*$", "", x_levels), ": ",
                                              sub("^.* \\| ", "", x_levels)), x_levels)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = FIG3_SMALL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.6),
          axis.text.y = element_text(size = 4.7),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5.2),
          legend.key.size = unit(7, "pt"))
  p_flip <- p_flip_bar / p_flip_dot +
    patchwork::plot_layout(heights = c(0.08, 1), guides = "keep")
  ggsave(file.path(FIG3_DIR, paste0(file_prefix, "_flipped.pdf")), p_flip,
         width = flip_width, height = flip_height, device = cairo_pdf)

  readr::write_csv(
    d %>% transmute(tissue = as.character(tissue), domain = as.character(domain),
                    marker_name, gene_id, marker_type, marker_class,
                    expected_domain = expected,
                    avg_expr = round(avg_expr, 4),
                    avg_scaled = round(avg_scaled, 4),
                    pct_exp = round(pct_exp, 4)),
    file.path(FIG3_DIR, paste0(file_prefix, "_values.csv")))
  invisible(list(plot = p_main, flipped = p_flip, values = d))
}

fig3a_main <- make_fig3a_marker_dotplot(
  fig3a_plot_df, FIG3A_TISSUES,
  "Fig3A_published_marker_dotplot_major_tissues",
  width = 5, height = 4, flip_width = 7, flip_height = 4.5)
fig3a_no_bud <- make_fig3a_marker_dotplot(
  fig3a_plot_df, setdiff(FIG3A_TISSUES, "Axillary bud"),
  "Fig3A_published_marker_dotplot_major_tissues_no_bud",
  width = 4.5, height = 4, flip_width = 6, flip_height = 4.5)

readr::write_csv(
  fig3_known_audit %>%
    filter(marker_name %in% fig3a_keep, tissue %in% FIG3A_TISSUES) %>%
    left_join(main_sel %>%
                select(gene_id, evidence_source, source_detail, ath_hit, func_anno),
              by = "gene_id") %>%
    transmute(marker_name, gene_id, marker_class, domain_accept = intended_domain,
              target_lab = vapply(intended_domain, fig3_short_target, character(1)),
              tissue, status, margin = round(margin, 4),
              best_target_domain = best_target, peak_domain,
              pct_detected_target = round(pct_target, 4), validated,
              evidence_source, source_detail,
              arabidopsis_best_hit = ath_hit, functional_annotation = func_anno),
  file.path(FIG3_DIR, "Fig3A_published_marker_dotplot_major_tissues_gene_list.csv"))

# ---------- 3.2 Fig. 3B: known markers in situ ----------
#
# Three markers spanning three concentric stem domains: KCS2 (epidermis),
# PP2-A10-2 (phloem), and 4CL (xylem). They are rendered on the stem
# cross-section series, where all three target domains are annotated and the
# spatial specificity claim can be tested directly.
#
# Open-ST capture is sparse, so a k=3 nearest-neighbour mean is applied FOR
# DISPLAY ONLY and every spot's raw value is written to the CSV alongside it.
# Smoothing lowers apparent contrast (see the supplementary comparison), so the
# statistics reported in the figure are conservative.

fig3_smooth_spatial <- function(coords, vals, k = 3) {
  if (nrow(coords) <= k) return(vals)
  nn <- FNN::get.knn(as.matrix(coords), k = k)$nn.index
  rowMeans(cbind(vals, matrix(vals[nn], nrow = length(vals))))
}

# Pull expression + harmonized domain + coordinates for one section set.
fig3_build_spatial <- function(obj, anno_col, sections, sel, k = 3) {
  dom_all <- harmonize_domain(as.character(obj@meta.data[[anno_col]]))
  names(dom_all) <- colnames(obj)
  purrr::map_dfr(sections, function(im) {
    co <- Seurat::GetTissueCoordinates(obj, image = im)
    co <- co[, c("x", "y")]
    cells <- rownames(co)
    purrr::pmap_dfr(list(sel$gene_id, sel$marker_name, sel$domain_lab),
      function(gi, mn, dl) {
        # SCT$data, the same normalized values the per-tissue caches and the DE
        # tables use -- reading Spatial$counts here would mix scales between the
        # maps and the statistics reported beside them
        e <- if (gi %in% rownames(obj[["SCT"]]$data))
               as.numeric(obj[["SCT"]]$data[gi, cells]) else rep(0, length(cells))
        tibble(x = co$x, y = co$y, cell = cells, sec = im,
               domain = dom_all[cells], raw = e,
               sm = fig3_smooth_spatial(co, e, k), marker = mn, dom_lab = dl)
      })
  })
}

# Section serials sit far apart in the capture coordinate frame, and each is
# offset diagonally from the last, so plotting them in native coordinates leaves
# most of the panel empty and shrinks the tissue to illegibility. Normalise each
# section to its own origin and tile left to right with a small gap, which is a
# change of layout only -- every per-spot value, domain assignment and statistic
# is computed before packing and is unaffected by it.
fig3_pack_sections <- function(df, order_secs, gap_frac = 0.06) {
  parts <- list(); xoff <- 0
  for (s in order_secs) {
    d <- df %>% filter(sec == s)
    if (!nrow(d)) next
    d$x <- d$x - min(d$x); d$y <- d$y - min(d$y)
    w <- max(d$x); d$x <- d$x + xoff
    parts[[s]] <- d
    xoff <- xoff + w * (1 + gap_frac)
  }
  bind_rows(parts)
}

FIG3B_SEL <- fig3a_known_sel %>%
  transmute(marker_name = as.character(marker_name),
            gene_id,
            marker_type = as.character(marker_type))

fig3b_genes <- unique(FIG3B_SEL$gene_id)
fig3b_sam_images <- "sam_A_s1"
fig3b_stem_images <- "stem_A_s1"
fig3b_bud_images <- 'bud_C_s2'
fig3b_petiole_cross_images <- "petiole_l4_cross1_s1"
  # "petiole_l4_cross1_s2"
fig3_spatial_theme <- theme(
  plot.title = element_text(size = 7, margin = margin(b = 4)),
  plot.margin = margin(8, 22, 8, 8),
  legend.title = element_text(size = 7),
  legend.text = element_text(size = 7),
  legend.key.height = unit(10, "pt"),
  legend.key.width = unit(6, "pt")
)

DefaultAssay(fig3_objs$sam) <- "SCT"
p_fig3b_sam <- SpatialFeaturePlot(
  fig3_objs$sam,
  features = intersect(fig3b_genes, rownames(fig3_objs$sam)),
  images = fig3b_sam_images,
  crop = TRUE,
  ncol = 2,
  pt.size.factor = 1.8,
  alpha = c(0.7, 0.9),
  min.cutoff = "q05",
  max.cutoff = "q95"
) & fig3_spatial_theme & theme(legend.position = "right")
save_pdf(p_fig3b_sam, file.path(FIG3_DIR, "Fig3B_known_marker_spatial_SAM_slice1.pdf"),
         width = 14, height = 10)

DefaultAssay(fig3_objs$stem_cross) <- "SCT"
p_fig3b_stem <- SpatialFeaturePlot(
  fig3_objs$stem_cross,
  features = intersect(fig3b_genes, rownames(fig3_objs$stem_cross)),
  images = fig3b_stem_images,
  crop = TRUE,
  ncol = 2,
  pt.size.factor = 5,
  alpha = c(0.7, 0.9),
  min.cutoff = "q05",
  max.cutoff = "q95"
) & fig3_spatial_theme & theme(legend.position = "right")
save_pdf(p_fig3b_stem, file.path(FIG3_DIR, "Fig3B_known_marker_spatial_Stem_slice1.pdf"),
         width = 14, height = 10)

DefaultAssay(fig3_objs$bud) <- "SCT"
p_fig3b_bud <- SpatialFeaturePlot(
  fig3_objs$bud,
  features = intersect(fig3b_genes, rownames(fig3_objs$bud)),
  images = fig3b_bud_images,
  crop = TRUE,
  ncol = 2,
  pt.size.factor = 4,
  alpha = c(0.7, 0.9),
  min.cutoff = "q05",
  max.cutoff = "q95"
) & fig3_spatial_theme & theme(legend.position = "right")
save_pdf(p_fig3b_bud, file.path(FIG3_DIR, "Fig3B_known_marker_spatial_Bud_all_slices.pdf"),
         width = 14, height = 12)

DefaultAssay(fig3_objs$petiole_cross) <- "SCT"
p_fig3b_petiole_cross <- SpatialFeaturePlot(
  fig3_objs$petiole_cross,
  features = intersect(fig3b_genes, rownames(fig3_objs$petiole_cross)),
  images = fig3b_petiole_cross_images,
  crop = TRUE,
  ncol = 2,
  pt.size.factor = 5,
  alpha = c(0.7, 0.9),
  min.cutoff = "q05",
  max.cutoff = "q95",
) & fig3_spatial_theme & theme(legend.position = "right")
save_pdf(p_fig3b_petiole_cross,
         file.path(FIG3_DIR, "Fig3B_known_marker_spatial_Petiole_cross_L4.pdf"),
         width = 14, height = 10)

fig3b_spatial_manifest <- tibble::tribble(
  ~output, ~tissue, ~object, ~images_used, ~n_markers_plotted, ~marker_gene_ids,
  "Fig3B_known_marker_spatial_SAM_slice1.pdf", "SAM", "sam",
  paste(fig3b_sam_images, collapse = ";"),
  length(intersect(fig3b_genes, rownames(fig3_objs$sam))),
  paste(intersect(fig3b_genes, rownames(fig3_objs$sam)), collapse = ";"),
  "Fig3B_known_marker_spatial_Stem_slice1.pdf", "Stem", "stem_cross",
  paste(fig3b_stem_images, collapse = ";"),
  length(intersect(fig3b_genes, rownames(fig3_objs$stem_cross))),
  paste(intersect(fig3b_genes, rownames(fig3_objs$stem_cross)), collapse = ";"),
  "Fig3B_known_marker_spatial_Bud_all_slices.pdf", "Axillary bud", "bud",
  paste(fig3b_bud_images, collapse = ";"),
  length(intersect(fig3b_genes, rownames(fig3_objs$bud))),
  paste(intersect(fig3b_genes, rownames(fig3_objs$bud)), collapse = ";"),
  "Fig3B_known_marker_spatial_Petiole_cross_L4.pdf", "Petiole cross", "petiole_cross",
  paste(fig3b_petiole_cross_images, collapse = ";"),
  length(intersect(fig3b_genes, rownames(fig3_objs$petiole_cross))),
  paste(intersect(fig3b_genes, rownames(fig3_objs$petiole_cross)), collapse = ";")
)
readr::write_csv(fig3b_spatial_manifest,
                 file.path(FIG3_DIR, "Fig3B_known_marker_spatial_manifest.csv"))

# Three serial sections of one stem, packed side by side for the de novo panel.
FIG3B_STEM_SECS <- c("stem_A_s1", "stem_A_s2", "stem_A_s3")

fig3_dom_pal <- setNames(
  scales::hue_pal(l = 55, c = 90)(length(FIG3_DOMAIN_LEVELS)), FIG3_DOMAIN_LEVELS)

fig3_map_domain <- function(d, title) {
  ggplot(d, aes(x, y, colour = domain)) +
    geom_point(size = 0.20, stroke = 0) +
    scale_colour_manual(values = fig3_dom_pal, na.value = "grey88", drop = TRUE) +
    coord_equal() + labs(title = title) + theme_void(base_size = FIG3_SMALL) +
    theme(legend.position = "none",
          plot.title = element_text(size = FIG3_SMALL, face = "bold", hjust = 0.5))
}

fig3_map_expr <- function(d, title) {
  ggplot(d, aes(x, y, colour = rel)) +
    geom_point(size = 0.20, stroke = 0) +
    scale_colour_gradientn(colours = c("grey90", "#FDD49E", "#FC8D59", "#B30000"),
                           limits = c(0, 1), name = "Rel. expr.") +
    coord_equal() + labs(title = title) + theme_void(base_size = FIG3_SMALL) +
    theme(legend.position = "none",
          plot.title = element_text(size = FIG3_SMALL, hjust = 0.5))
}

# Both keys are pulled off throwaway copies of the panels and placed once, so the
# maps themselves stay legend-free and every map keeps the same plotting area.
fig3_grab_legend <- function(p) cowplot::get_plot_component(
  p + theme(legend.position = "right", legend.text = element_text(size = 4.2),
            legend.title = element_text(size = 5),
            legend.key.width = unit(4, "pt"), legend.key.height = unit(8, "pt"),
            legend.key.size = unit(4, "pt")),
  "guide-box-right")

# Dotplot input straight from the per-tissue cache: domain-wise z-score of mean
# expression plus detection percentage. This helper is kept for supplementary
# panels where compact cached summaries are easier to lay out.
fig3_dotplot_from_cache <- function(cache, genes) {
  g <- intersect(genes, rownames(cache$avg))
  z <- t(scale(t(as.matrix(cache$avg[g, , drop = FALSE]))))
  z[!is.finite(z)] <- 0
  tidyr::expand_grid(gene_id = g, domain = colnames(cache$avg)) %>%
    mutate(avg_scaled = z[cbind(gene_id, domain)],
           pct_exp    = 100 * cache$pct[cbind(gene_id, domain)])
}

# ---------- 3.3 Fig. 3C: de novo marker dotplots ----------
#
# De novo markers are selected independently within each tissue/domain, excluding
# the curated known-marker set used in Fig. 3A. The all-tissue dotplots are
# screening panels; choose one anchor tissue for the main figure and move the
# remaining tissues to supplement.

fig3_denovo_pool <- function(de, cache, known_ids) {
  cl_col <- if ("cluster" %in% names(de)) "cluster" else "fig3_domain"
  if (!"gene_id" %in% names(de)) de$gene_id <- de$gene
  doms <- colnames(cache$avg)
  organelle_ids <- unique(c(read_feature_list("organelle_gene_ids.txt"),
                            read_feature_list("ChrPt_gene_ids.txt"),
                            read_feature_list("ChrMt_gene_ids.txt")))
  de %>%
    mutate(gene_symbol_for_filter = if ("gene" %in% names(.)) gene else gene_id) %>%
    filter(p_val_adj < 0.01, avg_log2FC > 0.5, pct.1 > 0.25,
           !gene_id %in% known_ids,
           !gene_id %in% organelle_ids,
           grepl("^PtXa", gene_id),
           !grepl("rrn|rRNA|trn|tRNA|gene-rrn|gene-trn|cds-", gene_id, ignore.case = TRUE),
           !grepl("rrn|rRNA|trn|tRNA|gene-rrn|gene-trn|cds-", gene_symbol_for_filter, ignore.case = TRUE),
           .data[[cl_col]] %in% doms, gene_id %in% rownames(cache$avg)) %>%
    rowwise() %>%
    mutate(tgt_expr    = cache$avg[gene_id, as.character(.data[[cl_col]])],
           max_off     = max(cache$avg[gene_id, setdiff(doms, as.character(.data[[cl_col]]))]),
           max_off_dom = setdiff(doms, as.character(.data[[cl_col]]))[
                           which.max(cache$avg[gene_id, setdiff(doms, as.character(.data[[cl_col]]))])],
           pct_t       = cache$pct[gene_id, as.character(.data[[cl_col]])]) %>%
    ungroup() %>%
    mutate(spec_log2fc = log2((tgt_expr + 0.05) / (max_off + 0.05)),
           passes_specificity = spec_log2fc > 0.5)
}

FIG3C_KEYS <- c("sam", "bud", "stem_cross", "petiole_cross")
FIG3C_TOP_N_PER_DOMAIN <- 4

fig3c_denovo <- lapply(FIG3C_KEYS, function(k) {
  de <- fig3_run_de(fig3_objs[[k]], k, FIG3_OBJ_SPEC[[k]]$anno)
  pool <- fig3_denovo_pool(de, fig3_caches[[k]], main_sel$gene_id)
  cl <- if ("cluster" %in% names(pool)) "cluster" else "fig3_domain"
  top <- pool %>%
    group_by(.data[[cl]]) %>%
    filter(if (any(passes_specificity)) passes_specificity else TRUE) %>%
    arrange(desc(passes_specificity), desc(pct_t), desc(tgt_expr),
            desc(spec_log2fc), desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = FIG3C_TOP_N_PER_DOMAIN) %>%
    mutate(selection_reason = ifelse(passes_specificity,
                                     "specific",
                                     "best_available_high_expression")) %>%
    ungroup()
  doms <- intersect(FIG3_DOMAIN_LEVELS, colnames(fig3_caches[[k]]$avg))
  ord <- top %>%
    mutate(dom = factor(as.character(.data[[cl]]), levels = doms)) %>%
    arrange(dom, desc(passes_specificity), desc(pct_t), desc(tgt_expr),
            desc(spec_log2fc), desc(avg_log2FC)) %>%
    pull(gene_id) %>%
    unique()
  list(tissue = FIG3_OBJ_SPEC[[k]]$tissue, obj_key = k, pool = pool,
       top = top, cl_col = cl, domains = doms, ord = ord,
       n_spots = sum(fig3_caches[[k]]$n_spots))
})
names(fig3c_denovo) <- FIG3C_KEYS

fig3c_pool <- purrr::map_dfr(fig3c_denovo, function(s)
  s$pool %>% mutate(tissue = s$tissue, obj_key = s$obj_key,
                    cluster = as.character(.data[[s$cl_col]])))
fig3c_top <- purrr::map_dfr(fig3c_denovo, function(s)
  s$top %>% mutate(tissue = s$tissue, obj_key = s$obj_key,
                   cluster = as.character(.data[[s$cl_col]])))
fig3c_dom_ord <- intersect(FIG3_DOMAIN_LEVELS, unique(as.character(fig3c_top$cluster)))

fig3c_dotplot_manifest <- purrr::imap_dfr(fig3c_denovo, function(s, k) {
  genes <- intersect(s$ord, rownames(fig3_caches[[k]]$avg))
  d <- fig3_dotplot_from_cache(fig3_caches[[k]], genes) %>%
    filter(gene_id %in% genes) %>%
    mutate(gene_id = factor(gene_id, levels = rev(genes)),
           domain = factor(as.character(domain), levels = s$domains),
           pct_plot = pmin(pct_exp, 75))
  p_vertical <- ggplot(d, aes(domain, gene_id)) +
    geom_point(aes(size = pct_plot, colour = avg_scaled), shape = 16, stroke = 0, alpha = 0.95) +
    scale_colour_gradient2(low = "#3B4CC0", mid = "grey92", high = "#B40426",
                           midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                           name = "Scaled\nexpr.") +
    scale_size_area(max_size = 4.8, limits = c(0, 75),
                    breaks = c(0, 25, 50, 75), name = "% spots") +
    labs(x = NULL, y = NULL, title = s$tissue) +
    theme_classic(base_size = FIG3_SMALL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
          axis.text.y = element_text(size = 5.8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.key.size = unit(8, "pt"),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5.2),
          plot.title = element_text(size = 8, face = "bold"))
  vertical_file <- file.path(
    FIG3_DIR,
    paste0("Fig3C_de_novo_marker_dotplot_", k, "_vertical.pdf")
  )
  save_pdf(p_vertical, vertical_file, width = 7.5,
           height = max(5, 0.19 * length(genes) + 1.8))

  tibble(
    tissue = s$tissue,
    object = k,
    vertical_output = basename(vertical_file),
    n_markers_plotted = length(genes),
    marker_gene_ids = paste(genes, collapse = ";")
  )
})
readr::write_csv(fig3c_dotplot_manifest,
                 file.path(FIG3_DIR, "Fig3C_de_novo_marker_dotplot_manifest.csv"))

readr::write_csv(
  fig3c_top %>% transmute(tissue, obj_key, gene_id,
                          domain = as.character(cluster),
                          avg_log2FC = round(avg_log2FC, 3), p_val_adj,
                          pct.1 = round(pct.1, 3), pct.2 = round(pct.2, 3),
                          spec_log2fc = round(spec_log2fc, 3),
                          target_mean = round(tgt_expr, 4),
                          best_off_mean = round(max_off, 4),
                          best_off_domain = max_off_dom,
                          pct_detected_target = round(pct_t, 3),
                          passes_specificity,
                          selection_reason,
                          in_known_marker_set = gene_id %in% main_sel$gene_id),
  file.path(FIG3_DIR, "Fig3C_de_novo_marker_gene_list_all_tissues.csv"))
readr::write_csv(fig3c_pool %>% mutate(across(where(is.numeric), ~round(., 4))),
                 file.path(FIG3_DIR, "Fig3C_de_novo_all_passing_candidates_all_tissues.csv"))

# ---------- 3.4 Fig. 3D: de novo markers in situ ----------

fig3d_sam_genes <- intersect(fig3c_denovo$sam$ord, rownames(fig3_objs$sam))
fig3d_stem_genes <- intersect(fig3c_denovo$stem_cross$ord, rownames(fig3_objs$stem_cross))
fig3d_bud_genes <- intersect(fig3c_denovo$bud$ord, rownames(fig3_objs$bud))
fig3d_petiole_cross_genes <- intersect(fig3c_denovo$petiole_cross$ord, rownames(fig3_objs$petiole_cross))

DefaultAssay(fig3_objs$sam) <- "SCT"
p_fig3d_sam <- SpatialFeaturePlot(
  fig3_objs$sam, features = fig3d_sam_genes, images = fig3b_sam_images,
  crop = TRUE, ncol = 2, pt.size.factor = 1.6,
  min.cutoff = "q05", max.cutoff = "q95"
) & fig3_spatial_theme & theme(legend.position = "right")
save_pdf(p_fig3d_sam, file.path(FIG3_DIR, "Fig3D_de_novo_marker_spatial_SAM_slice1.pdf"),
         width = 10, height = 17)

DefaultAssay(fig3_objs$stem_cross) <- "SCT"
p_fig3d_stem <- SpatialFeaturePlot(
  fig3_objs$stem_cross, features = fig3d_stem_genes, images = fig3b_stem_images,
  crop = TRUE, ncol = 2, pt.size.factor = 5,
  min.cutoff = "q05", max.cutoff = "q95"
) & fig3_spatial_theme & theme(legend.position = "right")
save_pdf(p_fig3d_stem, file.path(FIG3_DIR, "Fig3D_de_novo_marker_spatial_Stem_slice1.pdf"),
         width = 10, height = 17)

DefaultAssay(fig3_objs$bud) <- "SCT"
p_fig3d_bud <- SpatialFeaturePlot(
  fig3_objs$bud, features = fig3d_bud_genes, images = fig3b_bud_images,
  crop = TRUE, ncol = 2, pt.size.factor = 3,
  min.cutoff = "q05", max.cutoff = "q95"
) & fig3_spatial_theme & theme(legend.position = "right")
save_pdf(p_fig3d_bud, file.path(FIG3_DIR, "Fig3D_de_novo_marker_spatial_Bud.pdf"),
         width = 10, height = 17)

DefaultAssay(fig3_objs$petiole_cross) <- "SCT"
p_fig3d_petiole_cross <- SpatialFeaturePlot(
  fig3_objs$petiole_cross, features = fig3d_petiole_cross_genes,
  images = fig3b_petiole_cross_images,
  crop = TRUE, ncol = 2, pt.size.factor = 5,
  min.cutoff = "q05", max.cutoff = "q95"
) & fig3_spatial_theme & theme(legend.position = "right")
save_pdf(p_fig3d_petiole_cross,
         file.path(FIG3_DIR, "Fig3D_de_novo_marker_spatial_Petiole_cross.pdf"),
         width = 10, height = 17)

fig3d_pick <- fig3c_top %>%
  group_by(tissue, cluster) %>%
  slice_max(spec_log2fc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(tissue, marker_name = fig3_short_id(gene_id), gene_id,
            domain_lab = as.character(cluster), spec_log2fc)
fig3d_stats <- fig3d_pick %>%
  left_join(fig3c_top %>%
              transmute(tissue, gene_id, marker = fig3_short_id(gene_id),
                        dom_lab = as.character(cluster),
                        ratio = (tgt_expr + 0.05) / (max_off + 0.05),
                        ratio_min20 = ratio,
                        frac_detected = pct_t),
            by = c("tissue", "gene_id", "marker_name" = "marker",
                   "domain_lab" = "dom_lab"))
readr::write_csv(
  tibble::tribble(
    ~output, ~tissue, ~object, ~images_used, ~n_markers_plotted, ~marker_gene_ids,
    "Fig3D_de_novo_marker_spatial_SAM_slice1.pdf", "SAM", "sam",
    paste(fig3b_sam_images, collapse = ";"), length(fig3d_sam_genes),
    paste(fig3d_sam_genes, collapse = ";"),
    "Fig3D_de_novo_marker_spatial_Stem_slice1.pdf", "Stem", "stem_cross",
    paste(fig3b_stem_images, collapse = ";"), length(fig3d_stem_genes),
    paste(fig3d_stem_genes, collapse = ";"),
    "Fig3D_de_novo_marker_spatial_Bud.pdf", "Axillary bud", "bud",
    paste(fig3b_bud_images, collapse = ";"), length(fig3d_bud_genes),
    paste(fig3d_bud_genes, collapse = ";"),
    "Fig3D_de_novo_marker_spatial_Petiole_cross.pdf", "Petiole cross", "petiole_cross",
    paste(fig3b_petiole_cross_images, collapse = ";"), length(fig3d_petiole_cross_genes),
    paste(fig3d_petiole_cross_genes, collapse = ";")
  ),
  file.path(FIG3_DIR, "Fig3D_de_novo_marker_spatial_manifest.csv"))
readr::write_csv(fig3d_stats %>% mutate(across(where(is.numeric), ~round(., 4))),
                 file.path(FIG3_DIR, "Fig3D_spatial_map_statistics.csv"))

# ---------- 3.5 Supplementary marker atlas ----------
#
# Four supplements, all keyed to the main panels:
#   S1  per-tissue facets of the FULL 29-marker audit (superset of 3A)
#   S2  de novo dotplots for the four tissues not shown in 3C
#   S3  de novo yield per tissue - how many markers each domain supports
#   S4  unsmoothed spatial maps + smoothed/unsmoothed contrast comparison

# S1: same encoding as 3A but every curated marker, including the three that are
# testable in only one tissue.
fig3s_facets <- fig3_known_audit %>%
  mutate(class_s = factor(gsub(" / ", "/\n", marker_class),
                          levels = gsub(" / ", "/\n", FIG3A_CLASS_ORDER)),
         row_lab = paste0(marker_name, "  (",
                          vapply(intended_domain, fig3_short_target, character(1)), ")"),
         tissue_s = factor(tissue, levels = c("SAM", "Axillary bud", "Stem",
                                              "Petiole cross", "Petiole longitudinal")),
         margin_c = pmax(pmin(margin, 2), -2),
         state = case_when(status == "domain_not_annotated" ~ "Domain not annotated",
                           status == "gene_absent"          ~ "Gene not detected",
                           TRUE                             ~ NA_character_)) %>%
  arrange(class_s, desc(marker_name))
fig3s_facets$row_lab <- factor(fig3s_facets$row_lab, levels = unique(fig3s_facets$row_lab))

p_figs1 <- ggplot(fig3s_facets %>% filter(status == "tested"),
                  aes(tissue_s, row_lab, fill = margin_c)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#3B4CC0", mid = "grey96", high = "#B40426", midpoint = 0,
    limits = c(-2, 2), oob = scales::squish, name = "Specificity\nmargin (log2)",
    guide = guide_colourbar(barwidth = unit(4, "pt"), barheight = unit(26, "pt"),
                            title.position = "top")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = fig3s_facets %>% filter(!is.na(state)),
            aes(fill = state), colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Domain not annotated" = "grey80",
                               "Gene not detected"    = "grey55"), name = NULL,
                    guide = guide_legend(keywidth = unit(6, "pt"),
                                         keyheight = unit(6, "pt"))) +
  geom_point(data = fig3s_facets %>% filter(validated),
             aes(tissue_s, row_lab), colour = "white", size = 0.5, inherit.aes = FALSE) +
  facet_grid(class_s ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) + theme_bw(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.4),
        axis.text.y = element_text(size = 5),
        strip.text.y.left = element_text(angle = 0, size = 4.8, face = "bold",
                                         lineheight = 0.9, margin = margin(1, 2, 1, 2)),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey94", colour = NA),
        panel.grid = element_blank(), panel.spacing = unit(1.5, "pt"),
        legend.text = element_text(size = 4.8), legend.title = element_text(size = 5),
        legend.box.spacing = unit(3, "pt"))

# facet strips default to the width of the longest class label, which wastes
# roughly a third of the canvas; fix them at a set width on the assembled grob
figs1_grob <- ggplotGrob(p_figs1)
figs1_grob$widths[unique(figs1_grob$layout$l[grep("^strip-l", figs1_grob$layout$name)])] <-
  unit(34, "pt")
ggsave(file.path(FIG3_SUPP, "FigS_known_marker_per_tissue_facets.pdf"),
       figs1_grob, width = 3.9, height = 5.0, device = cairo_pdf)

readr::write_csv(
  fig3_known_audit %>%
    left_join(main_sel %>%
                select(gene_id, evidence_source, source_detail, ath_hit, func_anno),
              by = "gene_id") %>%
    transmute(marker_name, gene_id, marker_class, intended_domain, tissue, status,
              margin_log2 = round(margin, 3), best_target_domain = best_target,
              peak_domain, pct_detected_target = round(pct_target, 3), validated,
              plotted_in_Fig3A = marker_name %in% fig3a_keep,
              evidence_source, source_detail,
              arabidopsis_best_hit = ath_hit, functional_annotation = func_anno),
  file.path(FIG3_SUPP, "FigS_known_marker_full_audit_all_tissues.csv"))

# S2 + S3: the four remaining tissues. Same pool/rank logic as 3C so the panels
# are directly comparable; DE tables are cached per tissue.
FIG3S_KEYS <- c("sam", "bud", "petiole_cross", "petiole_long")

fig3s_denovo <- lapply(FIG3S_KEYS, function(k) {
  de   <- fig3_run_de(fig3_objs[[k]], k, FIG3_OBJ_SPEC[[k]]$anno)
  pool <- fig3_denovo_pool(de, fig3_caches[[k]], main_sel$gene_id)
  cl   <- if ("cluster" %in% names(pool)) "cluster" else "fig3_domain"
  top  <- pool %>% group_by(.data[[cl]]) %>%
    slice_max(spec_log2fc, n = 4, with_ties = FALSE) %>% ungroup()
  doms <- intersect(FIG3_DOMAIN_LEVELS, colnames(fig3_caches[[k]]$avg))
  ord  <- top %>% mutate(dom = factor(as.character(.data[[cl]]), levels = doms)) %>%
    arrange(dom, desc(spec_log2fc)) %>% pull(gene_id) %>% unique()
  list(tissue = FIG3_OBJ_SPEC[[k]]$tissue, pool = pool, top = top, cl_col = cl,
       domains = doms, n_spots = sum(fig3_caches[[k]]$n_spots),
       dp = fig3_dotplot_from_cache(fig3_caches[[k]], ord),
       ord = ord)
})
names(fig3s_denovo) <- FIG3S_KEYS

fig3s_dotplot <- function(s) {
  d <- s$dp %>% filter(gene_id %in% s$ord) %>%
    mutate(gene = factor(fig3_short_id(gene_id), levels = fig3_short_id(s$ord)),
           dom  = factor(as.character(domain), levels = s$domains))
  ggplot(d, aes(gene, dom)) +
    geom_point(aes(size = pct_exp, colour = avg_scaled)) +
    scale_size_area(max_size = 1.7, limits = c(0, 100), name = "% spots") +
    scale_colour_gradient2(low = "#3B4CC0", mid = "grey92", high = "#B40426",
                          midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                          name = "Scaled\nexpr.") +
    labs(x = NULL, y = NULL,
         title = sprintf("%s  (%d domains, %s spots)", s$tissue, length(s$domains),
                         format(s$n_spots, big.mark = ","))) +
    theme_bw(base_size = 5.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3.4),
          axis.text.y = element_text(size = 4.6),
          plot.title = element_text(size = 6, face = "bold"),
          panel.grid.major = element_line(linewidth = 0.15, colour = "grey92"),
          panel.grid.minor = element_blank(), legend.position = "none",
          plot.margin = margin(2, 2, 1, 2))
}

figs2_panels <- lapply(fig3s_denovo, fig3s_dotplot)
figs2_leg <- cowplot::get_plot_component(
  fig3s_dotplot(fig3s_denovo$sam) +
    theme(legend.position = "right", legend.text = element_text(size = 4.4),
          legend.title = element_text(size = 5), legend.key.size = unit(6, "pt")),
  "guide-box-right")
p_figs2 <- ((figs2_panels$sam / figs2_panels$bud /
             figs2_panels$petiole_cross / figs2_panels$petiole_long) |
            patchwork::wrap_elements(figs2_leg)) +
  patchwork::plot_layout(widths = c(1, 0.09))
ggsave(file.path(FIG3_SUPP, "FigS_de_novo_marker_dotplots_remaining_tissues.pdf"),
       p_figs2, width = 6.9, height = 6.2, device = cairo_pdf)

readr::write_csv(
  purrr::map_dfr(fig3s_denovo, function(s)
    s$top %>% transmute(tissue = s$tissue, gene_id, short_id = fig3_short_id(gene_id),
                        domain = as.character(.data[[s$cl_col]]),
                        avg_log2FC = round(avg_log2FC, 3), p_val_adj,
                        pct_target = round(pct.1, 3), pct_other = round(pct.2, 3),
                        specificity_log2FC = round(spec_log2fc, 3),
                        best_off_target_domain = max_off_dom,
                        in_known_marker_set = gene_id %in% main_sel$gene_id)),
  file.path(FIG3_SUPP, "FigS_de_novo_marker_gene_list_remaining_tissues.csv"))

# S3: yield table, including the anchor tissue from 3C
readr::write_csv(
  bind_rows(
    purrr::map_dfr(fig3s_denovo, function(s) {
      cl <- as.character(s$top[[s$cl_col]])
      tibble(tissue = s$tissue, n_domains = length(s$domains), n_spots = s$n_spots,
             n_de_novo_passing = nrow(s$pool), n_plotted = nrow(s$top),
             domains_with_marker = n_distinct(cl),
             domains_without_marker = paste(setdiff(s$domains, unique(cl)), collapse = "; "))
    }),
    tibble(tissue = FIG3_OBJ_SPEC$stem_cross$tissue,
           n_domains = length(fig3c_dom_ord),
           n_spots = sum(fig3_caches$stem_cross$n_spots),
           n_de_novo_passing = nrow(fig3c_pool), n_plotted = nrow(fig3c_top),
           domains_with_marker = n_distinct(as.character(fig3c_top$cluster)),
           domains_without_marker = paste(
             setdiff(intersect(FIG3_DOMAIN_LEVELS, colnames(fig3_caches$stem_cross$avg)),
                     unique(as.character(fig3c_top$cluster))), collapse = "; "))) %>%
    arrange(match(tissue, c("SAM", "Axillary bud", "Stem",
                            "Petiole cross", "Petiole longitudinal"))) %>%
    mutate(criteria = paste("padj<0.01, avg_log2FC>0.5, pct.1>0.25,",
                            "specificity log2FC>0.5 vs best off-target domain;",
                            "top up to 4 per domain plotted")),
  file.path(FIG3_SUPP, "FigS_de_novo_marker_yield_all_tissues.csv"))

# S4: de novo spatial maps now use Seurat SpatialFeaturePlot directly, without
# coordinate repacking or kNN smoothing. Keep a small audit table for traceability.
readr::write_csv(
  fig3d_stats %>%
    transmute(tissue, marker = marker_name, gene_id, claimed_domain = domain_lab,
              target_vs_best_off_ratio = round(ratio, 3),
              frac_detected = round(frac_detected, 3),
              note = "Spatial maps are direct Seurat SpatialFeaturePlot outputs; no kNN smoothing or coordinate repacking."),
  file.path(FIG3_SUPP, "FigS_de_novo_spatial_smoothing_comparison.csv"))

# ---------- 3.6 QC across every plotted gene ----------
#
# One row per gene that appears anywhere in Fig. 3, with its detection profile in
# all five tissues. This is the table to check before submission: it catches
# genes that are plotted but barely detected, and genes whose peak domain is
# inconsistent between tissues.

fig3_panel_map <- bind_rows(
  fig3_known_audit %>% filter(marker_name %in% fig3a_keep) %>%
    distinct(gene_id, marker_name) %>% mutate(panel = "Fig3A"),
  FIG3B_SEL  %>% distinct(gene_id, marker_name) %>% mutate(panel = "Fig3B"),
  fig3c_top  %>% distinct(gene_id) %>% mutate(marker_name = NA_character_, panel = "Fig3C"),
  fig3d_pick %>% distinct(gene_id, marker_name) %>% mutate(panel = "Fig3D"))

fig3_gene_profile <- purrr::map_dfr(unique(fig3_panel_map$gene_id), function(gi) {
  purrr::map_dfr(names(FIG3_OBJ_SPEC), function(k) {
    cc <- fig3_caches[[k]]
    if (!gi %in% rownames(cc$avg))
      return(tibble(gene_id = gi, tissue = FIG3_OBJ_SPEC[[k]]$tissue, in_matrix = FALSE,
                    peak_domain = NA_character_, max_pct = NA_real_))
    v <- cc$avg[gi, ]
    tibble(gene_id = gi, tissue = FIG3_OBJ_SPEC[[k]]$tissue, in_matrix = TRUE,
           peak_domain = colnames(cc$avg)[which.max(v)], max_pct = max(cc$pct[gi, ]))
  })
})

fig3_qc <- fig3_panel_map %>%
  group_by(gene_id) %>%
  summarise(panels = paste(sort(unique(panel)), collapse = ";"),
            label  = dplyr::first(na.omit(marker_name)), .groups = "drop") %>%
  mutate(gene_class = ifelse(gene_id %in% main_sel$gene_id,
                             "known (literature-curated)", "de novo (data-derived)")) %>%
  left_join(fig3_gene_profile %>% group_by(gene_id) %>%
              summarise(n_tissues_in_matrix = sum(in_matrix),
                        n_tissues_detected  = sum(in_matrix & max_pct >= 0.05, na.rm = TRUE),
                        max_pct_any_tissue  = max(max_pct, na.rm = TRUE),
                        peak_domains = paste(unique(na.omit(peak_domain)), collapse = ";"),
                        .groups = "drop"), by = "gene_id") %>%
  left_join(main_sel %>%
              select(gene_id, marker_class, intended_domain = domain_accept,
                     evidence_source, source_detail, ath_hit, func_anno), by = "gene_id") %>%
  left_join(fig3c_top %>% select(gene_id, denovo_domain = cluster,
                                 denovo_spec_log2fc = spec_log2fc,
                                 denovo_padj = p_val_adj, denovo_pct_target = pct_t),
            by = "gene_id") %>%
  left_join(fig3_known_audit %>% filter(status == "tested") %>% group_by(gene_id) %>%
              summarise(n_tissues_testable = n(),
                        n_tissues_validated = sum(validated), .groups = "drop"),
            by = "gene_id")

readr::write_csv(
  fig3_qc %>% transmute(
    gene_id, gene_label = coalesce(label, fig3_short_id(gene_id)),
    gene_class, panels, marker_class, intended_domain,
    n_tissues_in_expression_matrix = n_tissues_in_matrix,
    n_tissues_detected_ge5pct = n_tissues_detected,
    max_pct_detected_any_tissue = round(max_pct_any_tissue, 3),
    peak_domain_per_tissue = peak_domains,
    n_tissues_testable, n_tissues_validated,
    denovo_domain = as.character(denovo_domain),
    denovo_spec_log2fc = round(denovo_spec_log2fc, 3),
    denovo_padj, denovo_pct_target = round(denovo_pct_target, 3),
    denovo_in_known_set = gene_id %in% main_sel$gene_id,
    evidence_source, source_detail,
    arabidopsis_best_hit = ath_hit, functional_annotation = func_anno),
  file.path(FIG3_DIR, "Fig3_QC_all_plotted_genes.csv"))

fig3_all_gt1 <- function(x) length(x) > 0 && all(!is.na(x) & x > 1)
fig3_qc_checks <- tibble(
  check = c(
    "all curated known markers validate in at least one testable tissue",
    "every Fig3A marker is testable in at least two tissues",
    "known-marker Seurat spatial maps were generated for requested tissues",
    "de novo spatial target-domain contrast is >1",
    "de novo panels do not reuse curated known-marker genes",
    "every plotted gene is detected in at least one tissue"
  ),
  passed = c(
    all(fig3_known_audit %>% group_by(marker_name) %>%
          summarise(v = sum(validated), .groups = "drop") %>% pull(v) >= 1),
    all(fig3_known_audit %>% filter(marker_name %in% fig3a_keep) %>%
          group_by(marker_name) %>% summarise(n = sum(status == "tested"),
                                              .groups = "drop") %>% pull(n) >= 2),
    all(file.exists(file.path(FIG3_DIR, fig3b_spatial_manifest$output))) &&
      all(fig3b_spatial_manifest$n_markers_plotted > 0),
    fig3_all_gt1(c(fig3d_stats$ratio, fig3d_stats$ratio_min20)),
    length(intersect(fig3c_top$gene_id, main_sel$gene_id)) == 0,
    all(fig3_qc$n_tissues_in_matrix > 0) && all(fig3_qc$n_tissues_detected > 0)
  )
)
readr::write_csv(fig3_qc_checks, file.path(FIG3_DIR, "Fig3_QC_checks.csv"))
if (any(!fig3_qc_checks$passed)) {
  warning("Figure 3 QC checks not passed: ",
          paste(fig3_qc_checks$check[!fig3_qc_checks$passed], collapse = "; "))
}

writeLines(c(
  "# Figure 3 marker atlas design",
  "",
  "Core conclusion: the poplar spatial atlas recovers literature-supported cell-domain markers and nominates spatially resolved de novo markers.",
  "",
  "Panel A: published marker dotplot across SAM, axillary bud, stem, and petiole cross-section domains.",
  "Panel B: Seurat SpatialFeaturePlot contact sheets for all Fig. 3A published markers on SAM slice 1, stem slice 1, petiole longitudinal L4/L5, and petiole cross sections.",
  "Panel C: de novo stem marker dotplot using DotPlot_scCustom.",
  "Panel D: spatial maps for selected de novo stem markers on the same sections as panel B.",
  "",
  "Supplementary: full cross-tissue known-marker audit, de novo dotplots for other tissues, de novo marker-yield table, and smoothing audit.",
  "",
  "Suggested main-figure assembly: A and C as compact dotplots on the left/middle; B and D as spatial-map examples on the right or bottom. Keep the full cross-tissue conservation matrix in supplement."
), file.path(FIG3_DIR, "Fig3_design_and_legend.md"))

writeLines(c(
  "# Figure 3 QC report",
  "",
  paste("Known markers plotted in Fig3A:", paste(fig3a_keep, collapse = ", ")),
  paste("De novo stem markers plotted in Fig3C:", nrow(fig3c_top)),
  paste("De novo spatial markers plotted in Fig3D:", paste(fig3d_pick$marker_name, collapse = ", ")),
  "",
  "QC checks:",
  paste0("- ", fig3_qc_checks$check, ": ",
         ifelse(fig3_qc_checks$passed, "PASS", "REVIEW"))
), file.path(FIG3_DIR, "Fig3_QC_report.md"))

# MANUAL STITCH for Fig. 3:
#   Fig3A: published marker dotplot across major tissues.
#   Fig3B: manually selected examples from the Seurat SpatialFeaturePlot contact sheets.
#   Fig3C (4.2 x 2.0 in): de novo marker dotplot in stem.
#   Fig3D (5.2 x 2.1 in): de novo marker spatial examples.
#   Total budget ~7.6 x 5.2 in two-column layout. See
#   summary/Fig3_marker_atlas/Fig3_design_and_legend.md for the panel-size
#   budget and the draft legend.


# --------------------------------------- #
# S. Petiole developmental/spatial axis supplementary figure ####
# --------------------------------------- #

# Supplementary petiole axis: petiole cross and longitudinal annotation overview
plot_spatial_annotation(sobj_petiole_cross_ann, "Petiole cross annotation", "summary/Supplementary/FigS_petiole_axis/FigS_petiole_axis_A_cross_annotation.pdf")
plot_spatial_annotation(sobj_petiole_long_ann, "Petiole longitudinal annotation", "summary/Supplementary/FigS_petiole_axis/FigS_petiole_axis_B_longitudinal_annotation.pdf")

# Supplementary petiole axis: stage maps, if metadata exists
for (obj_nm in c("sobj_petiole_cross_ann")) {
  obj <- get(obj_nm)
  if (!is.null(obj)) {
    stage_col <- get_meta_col(obj, c("stage", "Stage", "leaf_stage", "petiole_stage", "sample_stage"))
    if (!is.null(stage_col)) {
      save_pdf(DimPlot(obj, reduction = "umap", group.by = stage_col, label = TRUE) + ggtitle(paste(obj_nm, stage_col)),
               paste0("summary/Supplementary/FigS_petiole_axis/FigS_petiole_axis_C_", obj_nm, "_stage_UMAP.pdf"), width = 8, height = 6)
      save_pdf(SpatialDimPlot(obj, group.by = stage_col, crop = FALSE, ncol = 3) + plot_annotation(title = paste(obj_nm, stage_col)),
               paste0("summary/Supplementary/FigS_petiole_axis/FigS_petiole_axis_D_", obj_nm, "_stage_spatial.pdf"), width = 14, height = 8)
    } else {
      message("No stage metadata found for ", obj_nm, ". Manual panel may use existing stage-specific figures.")
    }
  }
}

# Supplementary petiole axis: adaxial vs abaxial maps, if metadata exists
adab_col <- if (!is.null(sobj_petiole_cross_ann)) get_meta_col(sobj_petiole_cross_ann, c("adaxial_abaxial", "ad_ab", "region_adab", "manual_region")) else NULL
if (!is.null(adab_col)) {
  save_pdf(SpatialDimPlot(sobj_petiole_cross_ann, group.by = adab_col, crop = FALSE, ncol = 3),
           "summary/Supplementary/FigS_petiole_axis/FigS_petiole_axis_E_adaxial_abaxial_spatial.pdf", width = 14, height = 8)

  Idents(sobj_petiole_cross_ann) <- sobj_petiole_cross_ann[[adab_col]][, 1]
  adab_lvls <- levels(Idents(sobj_petiole_cross_ann))
  if (all(c("adaxial", "abaxial") %in% adab_lvls)) {
    petiole_adab_de <- FindMarkers(sobj_petiole_cross_ann, ident.1 = "adaxial", ident.2 = "abaxial",
                                   assay = "SCT", recorrect_umi = FALSE)
    petiole_adab_de$gene <- rownames(petiole_adab_de)
    write.csv(petiole_adab_de, "summary/Supplementary/FigS_petiole_axis/FigS_petiole_axis_F_adaxial_vs_abaxial_DE.csv", row.names = FALSE)
    top_adab_genes <- petiole_adab_de %>% filter(p_val_adj < 0.05) %>% arrange(desc(abs(avg_log2FC))) %>% slice_head(n = 12) %>% pull(gene)
    plot_spatial_features_safe(sobj_petiole_cross_ann, top_adab_genes,
                               "summary/Supplementary/FigS_petiole_axis/FigS_petiole_axis_G_adaxial_abaxial_top_genes.pdf")
  } else {
    message("Adaxial/abaxial labels are not exactly 'adaxial' and 'abaxial'. Edit Fig3D manually.")
  }
} else {
  message("No adaxial/abaxial metadata found. Use previous manually selected boundary plots if available.")
}

# MANUAL STITCH for Supplementary petiole axis:
#   Combine petiole cross/long annotations, stage maps, adaxial-abaxial maps,
#   and selected spatial DE genes. Use existing boundary-selection figures if
#   adaxial_abaxial metadata is not saved in the final object.


  },
  figure2 = function() {

local({
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(tidyr);library(tibble);library(ggplot2);library(patchwork);library(readr)})
workspace <- file.path(ATLAS_WORK_ROOT)
audit_dir <- file.path(ATLAS_CODE_ROOT, "resources", "marker_reference")
final_out <- file.path(workspace,'summary/Fig2_annotation/version2')
dir.create(final_out,recursive=TRUE,showWarnings=FALSE)
out <- tempfile('fig2_v2_stage_')
dir.create(out,recursive=TRUE)
source(file.path(ATLAS_CODE_ROOT, "R/figure_definitions.R"), local=TRUE)
stopifnot(exists('fig2c_reviewed_markers'))
original <- fig2c_reviewed_markers
markers <- original %>% mutate(change='Retained from provided figure script')
markers$marker_name[markers$marker_name=='4CL' & markers$gene_id=='PtXaAlbH.01G031600.v5.1'] <- '4CL1'
markers$change[markers$marker_name=='4CL1'] <- 'Label corrected; same locus'
# Preserve chr1 SEOR1 as its original locus, not as the published chr17 SEOR.
# Remove the unconfirmed chr3 PXY assignment from the main panel, but retain audit.
markers$plot_set[markers$gene_id=='PtXaAlbH.03G082400.v5.1'] <- 'supplement'
markers$change[markers$gene_id=='PtXaAlbH.03G082400.v5.1'] <- 'Unconfirmed PXY assignment; excluded from main panel'
add <- tribble(~tissue,~celltype,~marker_name,~gene_id,
'SAM','Xylem','4CL1','PtXaAlbH.01G031600.v5.1',
'Stem','Cambium','PXY','PtXaTreH.01G106100.v5.1',
'Stem','Xylem','4CL1','PtXaAlbH.01G031600.v5.1',
'Stem','Xylem','MAN4','PtXaTreH.06G095400.v5.1',
'Axillary bud','Xylem','MAN6','PtXaTreH.16G110600.v5.1',
'Axillary bud','Phloem','SEOR','PtXaAlbH.17G051000.v5.1',
'Petiole cross','Pith (Xylem)','4CL1','PtXaAlbH.01G031600.v5.1',
'Petiole cross','Pith (Xylem)','MAN6','PtXaTreH.16G110600.v5.1',
'Petiole longitudinal','Xylem1','4CL1','PtXaAlbH.01G031600.v5.1',
'Petiole longitudinal','Xylem1','MAN6','PtXaTreH.16G110600.v5.1',
'Petiole longitudinal','Phloem','SEOR','PtXaAlbH.17G051000.v5.1') %>% mutate(marker_role='Populus-supported marker',plot_set='main',notes='Selected using existing Populus locus and spatial expression audits',change='Added/promoted Populus-supported locus')
markers <- anti_join(markers,add,by=c('tissue','gene_id')) %>% bind_rows(add)
# Restore the complete marker set shown in the user's original stem panel.
restored <- tribble(~celltype,~marker_name,~gene_id,
'Epidermis','KCS2','PtXaAlbH.10G060600.v5.1',
'Epidermis','FDH','PtXaAlbH.06G180600.v5.1',
'Epidermis','LTL1','PtXaAlbH.19G028400.v5.1',
'Epidermis','LPTG1','PtXaAlbH.01G048200.v5.1',
'Cortex','PNP','PtXaAlbH.13G084300.v5.1',
'Cortex','PMEAMT','PtXaAlbH.15G031100.v5.1',
'Pith','HYDROLASE','PtXaAlbH.18G083900.v5.1',
'Pith','PSBXb','PtXaTreH.06G124900.v5.1',
'Pith','WRKY12','PtXaTreH.14G037200.v5.1',
'Cambium','ANT3','PtXaAlbH.14G005300.v5.1',
'Cambium','ACL5a','PtXaAlbH.06G184800.v5.1',
'Cambium','ANTL5-2','PtXaAlbH.18G067400.v5.1',
'Phloem','PIN3','PtXaTreH.10G093200.v5.1',
'Phloem','CLE41B','PtXaAlbH.02G198600.v5.1',
'Xylem','XCP1','PtXaAlbH.04G160300.v5.1',
'Xylem','XCP2','PtXaAlbH.05G203100.v5.1',
'Xylem','LAC4','PtXaAlbH.06G084500.v5.1',
'Xylem','4CL1','PtXaAlbH.01G031600.v5.1',
'Xylem','PHB','PtXaTreH.01G305500.v5.1') %>% mutate(tissue='Stem',plot_set='main',marker_role='Original figure domain marker',notes='Original figure nomenclature retained; 4CL label corrected to 4CL1. Domain enrichment is dataset-specific, not a claim of exclusive expression.',change='Retained/restored from original figure')
previous_stem <- markers %>% filter(tissue=='Stem',plot_set=='main')
markers <- bind_rows(filter(markers,tissue!='Stem'),restored,anti_join(previous_stem,restored,by='gene_id'))
stopifnot(all(restored$gene_id %in% markers$gene_id),all(previous_stem$gene_id %in% markers$gene_id))

# Compact main stem panel; retain the full candidate list in the selection audit.
stem_main_ids <- c(
 'PtXaAlbH.06G180600.v5.1','PtXaAlbH.10G060600.v5.1','PtXaAlbH.19G028400.v5.1',
 'PtXaAlbH.13G084300.v5.1','PtXaAlbH.15G031100.v5.1',
 'PtXaAlbH.18G083900.v5.1','PtXaTreH.06G124900.v5.1','PtXaTreH.14G037200.v5.1',
 'PtXaTreH.01G106100.v5.1','PtXaAlbH.18G067400.v5.1','PtXaAlbH.06G184800.v5.1',
 'PtXaAlbH.02G198600.v5.1','PtXaTreH.10G093200.v5.1','PtXaAlbH.06G071600.v5.1',
 'PtXaAlbH.01G031600.v5.1','PtXaAlbH.16G114900.v5.1','PtXaAlbH.05G203100.v5.1','PtXaAlbH.04G160300.v5.1')
if(!any(markers$tissue=='Stem' & markers$gene_id=='PtXaAlbH.16G114900.v5.1')) {
 markers <- bind_rows(markers,tibble(tissue='Stem',celltype='Xylem',marker_name='MAN6',gene_id='PtXaAlbH.16G114900.v5.1',marker_role='Populus-supported vessel candidate',plot_set='main',notes='Populus-supported locus with 17% stem xylem detection; paired with 4CL1 and the more widely detected XCP2.',change='Added to compact panel'))
}
markers <- markers %>% mutate(plot_set=ifelse(tissue=='Stem',ifelse(gene_id %in% stem_main_ids,'main','supplement'),plot_set))
stopifnot(sum(markers$tissue=='Stem' & markers$plot_set=='main')==18)
candidates <- read_csv(file.path(audit_dir,'poplar_marker_candidates.csv'),show_col_types=FALSE) %>% filter(marker %in% c('PXY','4CL1','XCP1','MAN4','MAN6','SEOR')) %>% mutate(gene_id=paste0(gene_id,'.v5.1'))
write_csv(markers,file.path(out,'Fig2_marker_selection_v2.csv'))
write_csv(candidates,file.path(out,'Populus_marker_evidence_v2.csv'))
files <- c('SAM'='sobj_sam_split_cleaned_res0.4_annotated_v2_2026.7.qs','Axillary bud'='sobj_bud_res0.5_final_v2_2026.7.24.qs','Stem'='sobj_stem_cross_AB_res0.6_final_v2_Seurat_RCTD_2026.7.26.qs','Petiole cross'='sobj_petiole_cross_res0.35_final_v1_2026.7.24.qs','Petiole longitudinal'='sobj_petiole_longitudinal_res0.15_final_v1_2026.7.24.qs')
save_plot <- function(p,stem,h,w=7.2) {
 ggsave(file.path(out,paste0(stem,'.pdf')),p,width=w,height=h,device=cairo_pdf)
 ggsave(file.path(out,paste0(stem,'.png')),p,width=w,height=h,dpi=220,bg='white')
 ggsave(file.path(out,paste0(stem,'.svg')),p,width=w,height=h,device=svglite::svglite)
}
all_data <- list(); counts <- list(); checks <- list()
for(tissue in names(files)) {
 message('Processing ',tissue)
 obj <- qread(file.path(workspace,'saved_obj',files[[tissue]]))
 if(tissue %in% c('SAM','Axillary bud')) {
  anno <- if(tissue=='SAM') sam_anno_v2 else bud_anno_v2
  labels <- setNames(anno$celltypes_v2,anno$cluster_label)[as.character(obj$seurat_clusters)]
  if(anyNA(labels)) stop('Unmapped clusters in ',tissue)
  obj$celltype <- unname(labels)
 } else obj$celltype <- as.character(obj$celltypes)
 stopifnot(!anyNA(obj$celltype))
 display <- if(tissue=='SAM') 'Shoot apex' else tissue
 slug <- gsub(' ','_',tolower(display))
 assay <- if('SCT' %in% Assays(obj)) 'SCT' else DefaultAssay(obj)
 DefaultAssay(obj) <- assay
 Idents(obj) <- 'celltype'
 m <- markers %>% filter(.data$tissue==!!tissue,plot_set=='main') %>% distinct(gene_id,.keep_all=TRUE)
 checks[[tissue]] <- m %>% mutate(present=gene_id %in% rownames(obj[[assay]]))
 if(any(!checks[[tissue]]$present)) stop('Main markers missing: ',tissue)
 # Original statistic: mean linearized normalized expression, row z-score across domains.
 features <- union(m$gene_id,candidates$gene_id)
 features <- intersect(features,rownames(obj[[assay]]))
 x <- GetAssayData(obj,assay=assay,layer='data')[features,,drop=FALSE]
 preferred <- switch(tissue,
  'SAM'=c('Shoot meristem','Proliferating','Leaf primordium','Epidermis','Mesophyll','Cortex','Pith','Xylem','Phloem'),
  'Stem'=c('Epidermis','Cortex','Pith','Cambium','Phloem','Xylem'),
  'Petiole cross'=c('Epidermis','Cortex','Inner cortex','Vasculature (Phloem)','Pith (Xylem)'),
  sort(unique(obj$celltype)))
 domains <- c(intersect(preferred,unique(obj$celltype)),setdiff(sort(unique(obj$celltype)),preferred))
 avg <- sapply(domains,function(g) Matrix::rowMeans(expm1(x[,obj$celltype==g,drop=FALSE])))
 pct <- sapply(domains,function(g) 100*Matrix::rowMeans(x[,obj$celltype==g,drop=FALSE]>0))
 z <- t(scale(t(avg))); z[!is.finite(z)] <- 0
 d <- expand_grid(gene_id=rownames(avg),celltype=colnames(avg)) %>% mutate(mean_expression=as.vector(t(avg)),pct_detected=as.vector(t(pct)),z_score=as.vector(t(z)),tissue=display,assay=assay)
 all_data[[tissue]] <- d
 counts[[tissue]] <- tibble(tissue=display,celltype=obj$celltype) %>% count(tissue,celltype,name='n_spots')
 draw_markers <- function(md,prefix) {
  md <- md %>% distinct(gene_id,.keep_all=TRUE) %>% filter(gene_id %in% features)
  # Sort rows by intended domain, preserving original within-domain ordering.
  md <- md %>% mutate(domain_order=match(tolower(celltype),tolower(domains))) %>% arrange(domain_order)
  labels <- setNames(paste(md$marker_name,md$gene_id,sep=' | '),md$gene_id)
  pd <- d %>% filter(gene_id %in% md$gene_id) %>% mutate(gene_id=factor(gene_id,levels=rev(md$gene_id)),celltype=factor(celltype,levels=domains))
  common <- theme_classic(base_size=9)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=8),axis.text.y=element_text(size=6.5,face='italic'),axis.title=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),plot.title=element_text(size=11))
  p <- ggplot(pd,aes(celltype,gene_id,fill=pmax(-2,pmin(2,z_score))))+geom_tile(color='white',linewidth=.25)+scale_y_discrete(labels=labels)+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-2,2),name='Scaled average\nexpression')+common+ggtitle(display)
  h <- max(4,1.7+nrow(md)*.19)
  save_plot(p,paste0(prefix,'_heatmap_v2'),h)
  p <- ggplot(pd,aes(celltype,gene_id,color=pmax(-2,pmin(2,z_score)),size=pct_detected))+geom_point()+scale_y_discrete(labels=labels)+scale_color_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-2,2),name='Scaled average\nexpression')+scale_size_area(max_size=4,limits=c(0,100),breaks=c(0,25,50,75,100),name='Spots detected (%)')+common+ggtitle(display)
  save_plot(p,paste0(prefix,'_dotplot_v2'),h)
 }
 draw_markers(m,paste0('Fig2C_',slug,'_poplar_markers'))
 cm <- candidates %>% transmute(gene_id,marker_name=marker,celltype=case_when(target_domain=='Xylem' & tissue=='Petiole cross'~'Pith (Xylem)',target_domain=='Xylem' & tissue=='Petiole longitudinal'~'Xylem1',TRUE~target_domain))
 if(tissue != 'Stem') draw_markers(cm,paste0('Supplement_',slug,'_all_Populus_candidates'))
 plot_spatial_annotation(obj,display,file.path(out,paste0('Fig2A_',slug,'_spatial_annotation_v2.pdf')),group_col='celltype',width=16,height=9)
 plot_umap_annotation(obj,display,file.path(out,paste0('Fig2B_',slug,'_UMAP_annotation_v2.pdf')),group_col='celltype')
 rm(obj,x);gc()
}
write_csv(bind_rows(all_data),file.path(out,'Fig2_marker_expression_source_data_v2.csv'))
write_csv(bind_rows(counts),file.path(out,'Fig2_domain_spot_counts_v2.csv'))
write_csv(bind_rows(checks),file.path(out,'Fig2_marker_presence_audit_v2.csv'))
writeLines(c('Figure 2 version 2: marker support for existing histology-informed tissue annotation.',
'Archetype: image plates plus quantitative marker grids; R-only rendering.',
'A/B: existing annotated objects and original plotting functions. No clustering or dimensional reduction was rerun; histology was not altered.',
'C: original main markers retained, with selected Populus-supported markers added; 4CL relabeled 4CL1 at the same locus. Unconfirmed chromosome-3 PXY removed from main panel.',
'Chromosome-1 SEOR1 retains its original name and locus; chromosome-17 SEOR is shown separately where supported.',
'MAN6 and XCP1 can have sparse detection. All Populus candidate loci are shown in supplementary heatmaps/dotplots across all tissues, including weaker loci.',
'Heatmap: SCT normalized data where available; expression linearized with expm1, averaged by annotated domain, then row z-scored across domains and clipped to [-2,2], matching original heatmap statistic.',
'Dot size: percentage of spots with nonzero expression in the same assay. Spot counts are descriptive, not independent biological replicates. No hypothesis tests or significance symbols.',
'Symbols and locus IDs italicized. PDF and SVG text remains editable. Individual panels require final assembly for submission.',
'Original outputs and saved objects are untouched. Source data and marker selection/evidence CSVs accompany these panels.'),file.path(out,'README_v2.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo_v2.txt'))

exports <- list.files(out,full.names=TRUE)
stopifnot(all(file.copy(exports,final_out,overwrite=TRUE)))
stopifnot(all(file.info(file.path(final_out,basename(exports)))$size == file.info(exports)$size))
message('Completed and verified ',length(exports),' files in ',final_out)

})

local({
suppressPackageStartupMessages({library(Seurat);library(qs);library(dplyr);library(tidyr);library(tibble);library(ggplot2);library(patchwork);library(readr)})
workspace <- file.path(ATLAS_WORK_ROOT)
audit_dir <- file.path(ATLAS_CODE_ROOT, "resources", "marker_reference")
final_out <- file.path(workspace,'summary/Fig2_annotation/version2')
dir.create(final_out,recursive=TRUE,showWarnings=FALSE)
out <- tempfile('fig2_v2_stage_')
dir.create(out,recursive=TRUE)
source(file.path(ATLAS_CODE_ROOT, "R/figure_definitions.R"), local=TRUE)
stopifnot(exists('fig2c_reviewed_markers'))
original <- fig2c_reviewed_markers
markers <- original %>% mutate(change='Retained from provided figure script')
markers$marker_name[markers$marker_name=='4CL' & markers$gene_id=='PtXaAlbH.01G031600.v5.1'] <- '4CL1'
markers$change[markers$marker_name=='4CL1'] <- 'Label corrected; same locus'
# Preserve chr1 SEOR1 as its original locus, not as the published chr17 SEOR.
# Remove the unconfirmed chr3 PXY assignment from the main panel, but retain audit.
markers$plot_set[markers$gene_id=='PtXaAlbH.03G082400.v5.1'] <- 'supplement'
markers$change[markers$gene_id=='PtXaAlbH.03G082400.v5.1'] <- 'Unconfirmed PXY assignment; excluded from main panel'
add <- tribble(~tissue,~celltype,~marker_name,~gene_id,
'SAM','Xylem','4CL1','PtXaAlbH.01G031600.v5.1',
'Stem','Cambium','PXY','PtXaTreH.01G106100.v5.1',
'Stem','Xylem','4CL1','PtXaAlbH.01G031600.v5.1',
'Stem','Xylem','MAN4','PtXaTreH.06G095400.v5.1',
'Axillary bud','Xylem','MAN6','PtXaTreH.16G110600.v5.1',
'Axillary bud','Phloem','SEOR','PtXaAlbH.17G051000.v5.1',
'Petiole cross','Pith (Xylem)','4CL1','PtXaAlbH.01G031600.v5.1',
'Petiole cross','Pith (Xylem)','MAN6','PtXaTreH.16G110600.v5.1',
'Petiole longitudinal','Xylem1','4CL1','PtXaAlbH.01G031600.v5.1',
'Petiole longitudinal','Xylem1','MAN6','PtXaTreH.16G110600.v5.1',
'Petiole longitudinal','Phloem','SEOR','PtXaAlbH.17G051000.v5.1') %>% mutate(marker_role='Populus-supported marker',plot_set='main',notes='Selected using existing Populus locus and spatial expression audits',change='Added/promoted Populus-supported locus')
markers <- anti_join(markers,add,by=c('tissue','gene_id')) %>% bind_rows(add)
# Restore the complete marker set shown in the user's original stem panel.
restored <- tribble(~celltype,~marker_name,~gene_id,
'Epidermis','KCS2','PtXaAlbH.10G060600.v5.1',
'Epidermis','FDH','PtXaAlbH.06G180600.v5.1',
'Epidermis','LTL1','PtXaAlbH.19G028400.v5.1',
'Epidermis','LPTG1','PtXaAlbH.01G048200.v5.1',
'Cortex','PNP','PtXaAlbH.13G084300.v5.1',
'Cortex','PMEAMT','PtXaAlbH.15G031100.v5.1',
'Pith','HYDROLASE','PtXaAlbH.18G083900.v5.1',
'Pith','PSBXb','PtXaTreH.06G124900.v5.1',
'Pith','WRKY12','PtXaTreH.14G037200.v5.1',
'Cambium','ANT3','PtXaAlbH.14G005300.v5.1',
'Cambium','ACL5a','PtXaAlbH.06G184800.v5.1',
'Cambium','ANTL5-2','PtXaAlbH.18G067400.v5.1',
'Phloem','PIN3','PtXaTreH.10G093200.v5.1',
'Phloem','CLE41B','PtXaAlbH.02G198600.v5.1',
'Xylem','XCP1','PtXaAlbH.04G160300.v5.1',
'Xylem','XCP2','PtXaAlbH.05G203100.v5.1',
'Xylem','LAC4','PtXaAlbH.06G084500.v5.1',
'Xylem','4CL1','PtXaAlbH.01G031600.v5.1',
'Xylem','PHB','PtXaTreH.01G305500.v5.1') %>% mutate(tissue='Stem',plot_set='main',marker_role='Original figure domain marker',notes='Original figure nomenclature retained; 4CL label corrected to 4CL1. Domain enrichment is dataset-specific, not a claim of exclusive expression.',change='Retained/restored from original figure')
previous_stem <- markers %>% filter(tissue=='Stem',plot_set=='main')
markers <- bind_rows(restored,anti_join(previous_stem,restored,by='gene_id'))
stopifnot(all(restored$gene_id %in% markers$gene_id),all(previous_stem$gene_id %in% markers$gene_id))

# Compact main stem panel; retain the full candidate list in the selection audit.
stem_main_ids <- c(
 'PtXaAlbH.06G180600.v5.1','PtXaAlbH.10G060600.v5.1','PtXaAlbH.19G028400.v5.1',
 'PtXaAlbH.13G084300.v5.1','PtXaAlbH.15G031100.v5.1',
 'PtXaAlbH.18G083900.v5.1','PtXaTreH.06G124900.v5.1','PtXaTreH.14G037200.v5.1',
 'PtXaTreH.01G106100.v5.1','PtXaAlbH.18G067400.v5.1','PtXaAlbH.06G184800.v5.1',
 'PtXaAlbH.02G198600.v5.1','PtXaTreH.10G093200.v5.1','PtXaAlbH.06G071600.v5.1',
 'PtXaAlbH.01G031600.v5.1','PtXaAlbH.16G114900.v5.1','PtXaAlbH.05G203100.v5.1','PtXaAlbH.04G160300.v5.1')
if(!any(markers$tissue=='Stem' & markers$gene_id=='PtXaAlbH.16G114900.v5.1')) {
 markers <- bind_rows(markers,tibble(tissue='Stem',celltype='Xylem',marker_name='MAN6',gene_id='PtXaAlbH.16G114900.v5.1',marker_role='Populus-supported vessel candidate',plot_set='main',notes='Populus-supported locus with 17% stem xylem detection; paired with 4CL1 and the more widely detected XCP2.',change='Added to compact panel'))
}
markers <- markers %>% mutate(plot_set=ifelse(tissue=='Stem',ifelse(gene_id %in% stem_main_ids,'main','supplement'),plot_set))
stopifnot(sum(markers$tissue=='Stem' & markers$plot_set=='main')==18)
candidates <- read_csv(file.path(audit_dir,'poplar_marker_candidates.csv'),show_col_types=FALSE) %>% filter(marker %in% c('PXY','4CL1','XCP1','MAN4','MAN6','SEOR')) %>% mutate(gene_id=paste0(gene_id,'.v5.1'))
write_csv(markers,file.path(out,'Fig2_marker_selection_v2.csv'))
write_csv(candidates,file.path(out,'Populus_marker_evidence_v2.csv'))
files <- c('SAM'='sobj_sam_split_cleaned_res0.4_annotated_v2_2026.7.qs','Axillary bud'='sobj_bud_res0.5_final_v2_2026.7.24.qs','Stem'='sobj_stem_cross_AB_res0.6_final_v2_Seurat_RCTD_2026.7.26.qs','Petiole cross'='sobj_petiole_cross_res0.35_final_v1_2026.7.24.qs','Petiole longitudinal'='sobj_petiole_longitudinal_res0.15_final_v1_2026.7.24.qs')
save_plot <- function(p,stem,h,w=7.2) {
 ggsave(file.path(out,paste0(stem,'.pdf')),p,width=w,height=h,device=cairo_pdf)
 ggsave(file.path(out,paste0(stem,'.png')),p,width=w,height=h,dpi=220,bg='white')
 ggsave(file.path(out,paste0(stem,'.svg')),p,width=w,height=h,device=svglite::svglite)
}
all_data <- list(); counts <- list(); checks <- list()
for(tissue in 'Stem') {
 message('Processing ',tissue)
 obj <- qread(file.path(workspace,'saved_obj',files[[tissue]]))
 if(tissue %in% c('SAM','Axillary bud')) {
  anno <- if(tissue=='SAM') sam_anno_v2 else bud_anno_v2
  labels <- setNames(anno$celltypes_v2,anno$cluster_label)[as.character(obj$seurat_clusters)]
  if(anyNA(labels)) stop('Unmapped clusters in ',tissue)
  obj$celltype <- unname(labels)
 } else obj$celltype <- as.character(obj$celltypes)
 stopifnot(!anyNA(obj$celltype))
 display <- if(tissue=='SAM') 'Shoot apex' else tissue
 slug <- gsub(' ','_',tolower(display))
 assay <- if('SCT' %in% Assays(obj)) 'SCT' else DefaultAssay(obj)
 DefaultAssay(obj) <- assay
 Idents(obj) <- 'celltype'
 m <- markers %>% filter(.data$tissue==!!tissue,plot_set=='main') %>% distinct(gene_id,.keep_all=TRUE)
 checks[[tissue]] <- m %>% mutate(present=gene_id %in% rownames(obj[[assay]]))
 if(any(!checks[[tissue]]$present)) stop('Main markers missing: ',tissue)
 # Original statistic: mean linearized normalized expression, row z-score across domains.
 features <- union(m$gene_id,candidates$gene_id)
 features <- intersect(features,rownames(obj[[assay]]))
 x <- GetAssayData(obj,assay=assay,layer='data')[features,,drop=FALSE]
 preferred <- switch(tissue,
  'SAM'=c('Shoot meristem','Proliferating','Leaf primordium','Epidermis','Mesophyll','Cortex','Pith','Xylem','Phloem'),
  'Stem'=c('Epidermis','Cortex','Pith','Cambium','Phloem','Xylem'),
  'Petiole cross'=c('Epidermis','Cortex','Inner cortex','Vasculature (Phloem)','Pith (Xylem)'),
  sort(unique(obj$celltype)))
 domains <- c(intersect(preferred,unique(obj$celltype)),setdiff(sort(unique(obj$celltype)),preferred))
 avg <- sapply(domains,function(g) Matrix::rowMeans(expm1(x[,obj$celltype==g,drop=FALSE])))
 pct <- sapply(domains,function(g) 100*Matrix::rowMeans(x[,obj$celltype==g,drop=FALSE]>0))
 z <- t(scale(t(avg))); z[!is.finite(z)] <- 0
 d <- expand_grid(gene_id=rownames(avg),celltype=colnames(avg)) %>% mutate(mean_expression=as.vector(t(avg)),pct_detected=as.vector(t(pct)),z_score=as.vector(t(z)),tissue=display,assay=assay)
 all_data[[tissue]] <- d
 counts[[tissue]] <- tibble(tissue=display,celltype=obj$celltype) %>% count(tissue,celltype,name='n_spots')
 draw_markers <- function(md,prefix) {
  md <- md %>% distinct(gene_id,.keep_all=TRUE) %>% filter(gene_id %in% features)
  # Sort rows by intended domain, preserving original within-domain ordering.
  md <- md %>% mutate(domain_order=match(tolower(celltype),tolower(domains))) %>% arrange(domain_order)
  labels <- setNames(paste(md$marker_name,md$gene_id,sep=' | '),md$gene_id)
  pd <- d %>% filter(gene_id %in% md$gene_id) %>% mutate(gene_id=factor(gene_id,levels=rev(md$gene_id)),celltype=factor(celltype,levels=domains))
  common <- theme_classic(base_size=9)+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=8),axis.text.y=element_text(size=6.5,face='italic'),axis.title=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),plot.title=element_text(size=11))
  p <- ggplot(pd,aes(celltype,gene_id,fill=pmax(-2,pmin(2,z_score))))+geom_tile(color='white',linewidth=.25)+scale_y_discrete(labels=labels)+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-2,2),name='Scaled average\nexpression')+common+ggtitle(display)
  h <- max(4,1.7+nrow(md)*.19)
  save_plot(p,paste0(prefix,'_heatmap_v2'),h)
  p <- ggplot(pd,aes(celltype,gene_id,color=pmax(-2,pmin(2,z_score)),size=pct_detected))+geom_point()+scale_y_discrete(labels=labels)+scale_color_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-2,2),name='Scaled average\nexpression')+scale_size_area(max_size=4,limits=c(0,100),breaks=c(0,25,50,75,100),name='Spots detected (%)')+common+ggtitle(display)
  save_plot(p,paste0(prefix,'_dotplot_v2'),h)
 }
 draw_markers(m,paste0('Fig2C_',slug,'_poplar_markers'))
 cm <- candidates %>% transmute(gene_id,marker_name=marker,celltype=case_when(target_domain=='Xylem' & tissue=='Petiole cross'~'Pith (Xylem)',target_domain=='Xylem' & tissue=='Petiole longitudinal'~'Xylem1',TRUE~target_domain))
 if(tissue != 'Stem') draw_markers(cm,paste0('Supplement_',slug,'_all_Populus_candidates'))
 rm(obj,x);gc()
}
write_csv(bind_rows(all_data),file.path(out,'Fig2_marker_expression_source_data_v2.csv'))
write_csv(bind_rows(counts),file.path(out,'Fig2_domain_spot_counts_v2.csv'))
write_csv(bind_rows(checks),file.path(out,'Fig2_marker_presence_audit_v2.csv'))
writeLines(c('Compact stem panel: 18 loci; epidermis 3, cortex 2, pith 3, cambium 3, phloem 3, xylem 4. XCP1 is included in the main panel. Unselected loci remain in the internal selection audit.',
'HYDROLASE and PSBXb are restored as empirical pith-enriched markers; literature support supplements rather than replaces informative dataset markers.',
'Original figure gene symbols retained except 4CL corrected to 4CL1. The original Alb FDH locus is retained; duplicate FDH alleles are omitted from the main panel.',
'Same saved stem object, SCT assay, six-domain order, linearized mean-expression calculation, row z-score and clipping at +/-2 as v2. No clustering was rerun.',
'The stem is presented in the main figure only; no supplementary stem marker panel is generated. Dot sizes show nonzero detection percentages.',
'PDF/SVG editable panels, PNG previews, source data and selected-marker audit accompany this script. Earlier outputs remain untouched.'),file.path(out,'README_v2.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo_v2.txt'))

# Keep other tissues in shared v2 source tables while replacing stem records.
for(nm in c('Fig2_marker_selection_v2.csv','Fig2_marker_expression_source_data_v2.csv','Fig2_domain_spot_counts_v2.csv','Fig2_marker_presence_audit_v2.csv')) {
 oldpath <- file.path(final_out,nm)
 if(file.exists(oldpath)) {
  old <- read_csv(oldpath,show_col_types=FALSE)
  updated <- read_csv(file.path(out,nm),show_col_types=FALSE)
  write_csv(bind_rows(filter(old,tissue!='Stem'),updated),file.path(out,nm))
 }
}
file.rename(file.path(out,'README_v2.txt'),file.path(out,'README_stem_update_v2.txt'))
file.rename(file.path(out,'R_sessionInfo_v2.txt'),file.path(out,'R_sessionInfo_stem_update_v2.txt'))
exports <- list.files(out,full.names=TRUE)
stopifnot(all(file.copy(exports,final_out,overwrite=TRUE)))
stopifnot(all(file.info(file.path(final_out,basename(exports)))$size == file.info(exports)$size))
message('Completed and verified ',length(exports),' files in ',final_out)

})

  },
  figure3 = function() {
# Cross-organ tissue-associated marker comparison; existing source values retained.
# Quantitative grid, four organ blocks, editable PDF/SVG and 300 dpi preview.
library(grid)
library(svglite)
library(ragg)
base <- file.path(ATLAS_WORK_ROOT, 'summary/Fig3_marker_atlas')
out <- file.path(base,'version2'); dir.create(out,showWarnings=FALSE)
d <- read.csv(file.path(base,'Fig3A_published_marker_dotplot_major_tissues_values.csv'),check.names=FALSE)
d <- d[!d$marker_name %in% c('PXY','AtHB8'),]
d$marker_name[d$marker_name=='ANTb'] <- 'ANT'
d$marker_name[d$marker_name=='4CL'] <- '4CL1'
d$marker_name[d$marker_name=='MAN6a/b'] <- 'MAN6'
d$domain[d$tissue=='Petiole cross' & d$domain=='Pith (Xylem)'] <- 'Xylem'
d$domain[d$tissue=='Petiole cross' & d$domain=='Vasculature (Phloem)'] <- 'Phloem'
d$tissue[d$tissue=='SAM'] <- 'Shoot apex'
organs <- c('Shoot apex','Axillary bud','Stem','Petiole cross')
genes <- c('KCS2','FDH','WRKY12','ATHB13','ANT','CYCB1;5','PP2-A10-2','SEOR1','4CL1','MAN6')
cols <- unique(d[c('tissue','domain')]); cols <- do.call(rbind,lapply(organs,function(t)cols[cols$tissue==t,])); cols$x <- seq_len(nrow(cols))
d$x <- cols$x[match(paste(d$tissue,d$domain),paste(cols$tissue,cols$domain))]
d$y <- 11-match(d$marker_name,genes)
stopifnot(nrow(d)==10*nrow(cols),!anyNA(d$x),!anyNA(d$y))
write.csv(d,file.path(out,'Fig3A_source_data_v2.csv'),row.names=FALSE)
# One drawing viewport provides identical boundaries for header strips and dots.
plotit <- function(){
 grid.newpage(); pushViewport(viewport(x=.025,y=.17,width=.95,height=.72,just=c('left','bottom'),xscale=c(-11,nrow(cols)+4),yscale=c(.5,11.7)))
 xx <- function(x)unit(x,'native'); yy <- function(y)unit(y,'native')
 tx <- function(s,x,y,size=10,...){grid.text(s,x=xx(x),y=yy(y),gp=gpar(fontfamily='Arial',fontsize=size),...)}
 pal <- colorRampPalette(c('#4D48B7','#EDEDED','#B60826'))(401)
 organ_colors <- c('#E467E9','#F7766D','#A3A500','#00ADEF')
 for(i in seq_along(organs)){
  cc <- cols$x[cols$tissue==organs[i]]; l <- min(cc)-.5; rr <- max(cc)+.5
  grid.rect(x=xx((l+rr)/2),y=yy(11.15),width=unit(rr-l,'native'),height=unit(.9,'native'),gp=gpar(fill=organ_colors[i],col=NA))
  tx(organs[i],(l+rr)/2,11.15,14)
  if(i>1)grid.lines(x=xx(c(l,l)),y=yy(c(.5,10.5)),gp=gpar(col='#AAAAAA',lwd=1))
 }
 groups <- c('Epidermis','Ground tissue','Meristem\nprimordium','Phloem','Xylem')
 for(i in 1:5){
  y <- 10.5-(i-.5)*2
  grid.rect(x=xx(-8.5),y=yy(y),width=unit(4.8,'native'),height=unit(1.98,'native'),gp=gpar(fill='#E0E0E0',col=NA))
  tx(groups[i],-8.5,y,12)
 }
 for(i in seq_along(genes)){
  id <- unique(d$gene_id[d$marker_name==genes[i]]);stopifnot(length(id)==1)
  y<-11-i
  grid.text(genes[i],x=xx(.4),y=yy(y+.13),just='right',gp=gpar(fontfamily='Arial',fontsize=10,fontface='italic'))
  grid.text(id,x=xx(.4),y=yy(y-.13),just='right',gp=gpar(fontfamily='Arial',fontsize=9,fontface='italic'))
 }
 for(i in seq_len(nrow(d))){
  col <- pal[round((max(-2,min(2,d$avg_scaled[i]))+2)*100)+1]
  grid.circle(x=xx(d$x[i]),y=yy(d$y[i]),r=unit(.65+2.15*sqrt(d$pct_exp[i]/100),'mm'),gp=gpar(fill=col,col=NA))
 }
 for(i in seq_len(nrow(cols)))tx(cols$domain[i],cols$x[i],.33,10,just='right',rot=48)
 lx<-nrow(cols)+2
 tx('Scaled\nexpression',lx,7.8,11)
 for(i in 1:200)grid.rect(x=xx(lx-.3),y=yy(5+(i-.5)/100),width=unit(.65,'native'),height=unit(.011,'native'),gp=gpar(fill=pal[2*i],col=NA))
 for(v in -2:2)tx(as.character(v),lx+.5,6+v/2,10)
 tx('% spots',lx,4.35,11)
 for(i in 1:5){v<-c(0,25,50,75,100)[i];y<-4-i*.48;grid.circle(x=xx(lx-.3),y=yy(y),r=unit(.65+2.15*sqrt(v/100),'mm'),gp=gpar(fill='#111111',col=NA));tx(v,lx+.55,y,10)}
 popViewport()
}
f<-file.path(out,'Fig3A_published_marker_dotplot_v2')
cairo_pdf(paste0(f,'.pdf'),width=14,height=10,family='Arial');plotit();dev.off()
svglite(paste0(f,'.svg'),width=14,height=10);plotit();dev.off()
agg_png(paste0(f,'.png'),width=14,height=10,units='in',res=300);plotit();dev.off()

  },
  figure4 = function() {

local({
suppressPackageStartupMessages({library(Seurat);library(qs);library(ggplot2);library(patchwork);library(dplyr);library(tidyr);library(monocle3)})
base <- file.path(ATLAS_WORK_ROOT)
out <- file.path(base,'summary/Fig4_trichome/version2');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x <- readRDS(atlas_input("saved_obj/figure4_plot_inputs.rds"))
list2env(x, environment())
s <- qread(file.path(base,'saved_obj/sobj_shoot_apex_trichome_scores_companion_2026-09-14.qs'))
export <- function(p,n,w=183,h=100){for(ext in c('pdf','svg','png'))ggsave(file.path(out,paste0(n,'.',ext)),p,width=w,height=h,units='mm',device=switch(ext,pdf=cairo_pdf,svg=svglite::svglite,png=ragg::agg_png),dpi=200,limitsize=FALSE,bg='white')}
clean <- function(p) p+labs(title=NULL,subtitle=NULL)+theme(plot.title=element_blank(),plot.subtitle=element_blank(),text=element_text(family='Arial',size=8))
pa <- clean(pA)+scale_x_discrete(labels=function(z){z[z=='14: Epidermis']<-'14: Epidermis (putative trichome initials)';z[z=='39: Epidermis']<-'39: Epidermis (putative developing trichomes)';z})
pB <- SpatialFeaturePlot(s,features=c('trichome_core_score','rctd_trichome_initial','rctd_trichome_developing'),images='sam_A_s1',crop=TRUE,combine=FALSE,pt.size.factor=1.75)
for(i in 1:3)pB[[i]]<-pB[[i]]+labs(title=c('Core trichome score','Predicted scRNA-seq cluster 14\nPutative trichome initials','Predicted scRNA-seq cluster 39\nPutative developing trichomes')[i],fill=c('Score','Predicted weight','Predicted weight')[i])+theme(legend.position='bottom',plot.title=element_text(size=8),legend.title=element_text(size=6),legend.text=element_text(size=6))+guides(fill=guide_colorbar(title.position='top',barwidth=grid::unit(20,'mm'),barheight=grid::unit(2,'mm')))
pb <- wrap_plots(pB,nrow=1)
# Existing inferred trajectory; no graph or pseudotime recomputation.
cds <- qread(file.path(base,'trajectory/SAM/monocle3/SAM_upward_leaf_branch_monocle3_cds.qs'))
u <- upward;rownames(u)<-u$spot
common <- intersect(colnames(cds),rownames(u)); stopifnot(length(common)>100); cds <- cds[,common]
colData(cds)$developmental_subcluster <- as.character(u[colnames(cds),'subcluster'])
pcgraph <- plot_cells(cds,color_cells_by='developmental_subcluster',label_cell_groups=FALSE,label_leaves=FALSE,label_roots=FALSE,label_branch_points=FALSE,show_trajectory_graph=TRUE,cell_size=.55)+scale_color_manual(values=subcluster_colors,name='Subcluster')+labs(colour='Subcluster',title='Developmental trajectory')+theme_classic(base_size=8)
s$developmental_subcluster <- factor(u[colnames(s),'subcluster'],levels=names(subcluster_colors))
pcmap <- SpatialDimPlot(s,group.by='developmental_subcluster',images='sam_A_s1',cols=subcluster_colors,crop=TRUE,pt.size.factor=1.75)+labs(title='Spatial subclusters',fill='Subcluster')+theme(legend.position='none',plot.title=element_text(size=9))
pc <- pcmap+pcgraph+plot_layout(widths=c(.7,1.3))
# Marker panel is predefined, not filtered by trajectory-test significance.
mk <- marker_tbl[marker_tbl$include_in_core,]
spots <- u$spot[order(u$monocle3_pseudotime_upward_leaf)]
bins <- dplyr::ntile(seq_along(spots),20)
a <- as.matrix(GetAssayData(s,assay='SCT',layer='data')[mk$spatial_feature,spots,drop=FALSE])
bm <- sapply(1:20,function(b)rowMeans(a[,bins==b,drop=FALSE]));z<-t(scale(t(bm)));z[!is.finite(z)]<-0
hm <- as.data.frame(as.table(z));names(hm)<-c('gene','bin','z');hm$bin<-rep(1:20,each=nrow(z));hm$label<-factor(paste0(mk$gene_name[match(hm$gene,mk$spatial_feature)],' | ',mk$gene_id[match(hm$gene,mk$spatial_feature)]),levels=rev(paste0(mk$gene_name,' | ',mk$gene_id)))
# Full trajectory-significant set; marker labels point to actual rows.
full_src <- file.path(out,'FigS_full_trajectory_gene_heatmap_source.csv')
if (TRUE) {
 dg <- read.csv(file.path(base,'tables/trajectory/SAM_upward_leaf_branch_dynamic_genes.csv'))
 ids <- intersect(dg$gene[!is.na(dg$q_value) & dg$q_value<.01 & dg$morans_I>0],rownames(s[['SCT']]))
 aa <- GetAssayData(s,assay='SCT',layer='data')[ids,spots,drop=FALSE]
 bb <- sapply(1:20,function(b)Matrix::rowMeans(aa[,bins==b,drop=FALSE]))
 zz <- t(scale(t(bb)));zz[!is.finite(zz)]<-0
 zz <- zz[order(max.col(zz,ties.method='first')),,drop=FALSE]
 hh <- data.frame(gene=rep(rownames(zz),20),row=rep(seq_len(nrow(zz)),20),bin=rep(1:20,each=nrow(zz)),z=as.vector(zz))
} else hh <- read.csv(full_src)
ann <- unique(hh[c('gene','row')]) |> inner_join(mk,by=c('gene'='spatial_feature')) |> arrange(row)
stopifnot(nrow(ann)==6)
ann$label <- paste0(ann$gene_name,' | ',ann$gene_id)
ann$label_y <- -ann$row
sep <- max(hh$row)*.065
for(i in seq_len(nrow(ann))[-1]) ann$label_y[i] <- min(ann$label_y[i],ann$label_y[i-1]-sep)
# Shift crowded labels upward if needed, preserving their order and leader lines.
if(min(ann$label_y) < -max(hh$row)) ann$label_y <- ann$label_y + (-max(hh$row)-min(ann$label_y))
pd <- ggplot(hh,aes(bin,-row,fill=z))+geom_raster()+
 geom_segment(data=ann,aes(x=20.6,xend=22,y=-row,yend=label_y),inherit.aes=FALSE,linewidth=.25,color='#555555')+
 geom_text(data=ann,aes(x=22.3,y=label_y,label=label),inherit.aes=FALSE,hjust=0,size=2.6,fontface='italic')+
 scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-2,2),oob=scales::squish)+
 scale_x_continuous(limits=c(.5,41),breaks=c(1,5,10,15,20),expand=c(0,0))+
 scale_y_continuous(expand=c(.005,0))+
 labs(x='Developmental pseudotime → (20 equal-count bins)',y=paste(max(hh$row),'trajectory-significant genes'),fill='Gene z-score')+
 theme_minimal(base_size=8)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.grid=element_blank(),legend.position='bottom')
write.csv(ann[c('gene','gene_name','row','label_y')],file.path(out,'Fig4d_highlighted_gene_rows_v2.csv'),row.names=FALSE)
write.csv(hh,file.path(out,'Fig4d_full_trajectory_heatmap_source_v2.csv'),row.names=FALSE)
write.csv(hm,file.path(out,'Fig4d_marker_pseudotime_source.csv'),row.names=FALSE)
dyn <- read.csv(file.path(base,'tables/trajectory/SAM_upward_leaf_branch_dynamic_genes.csv'))
write.csv(mk|>left_join(dyn,by=c('spatial_feature'='gene')),file.path(out,'Fig4d_marker_trajectory_test_audit.csv'),row.names=FALSE)
pe <- clean(pE)+labs(y='Core trichome score',x='Developmental pseudotime')
pf <- clean(pF)+scale_x_discrete(labels=function(z)sub('SAM-','Shoot apex-',z))+scale_y_discrete(labels=c('Core trichome score'='Core trichome score','Predicted scRNA-seq cluster 14 contribution (trichome initials)'='Predicted scRNA-seq cluster 14\n(putative trichome initials)','Predicted scRNA-seq cluster 39 contribution (developing trichomes)'='Predicted scRNA-seq cluster 39\n(putative developing trichomes)','Meristem-to-primordium pseudotime'='Developmental pseudotime'))
export(pa,'FigS_single_cell_reference_markers_v2',183,110);export(pb,'Fig4b_spatial_score_and_predictions_v2',183,85);export(pc,'Fig4a_subclusters_and_trajectory_v2',183,95);export(pd,'Fig4c_full_trajectory_heatmap_v2',183,190);export(pe,'Fig4d_core_score_pseudotime_v2',150,100);export(pf,'Fig4e_module_correlations_v2',210,85)
# Assemble with explicit tags so nested spatial panels do not acquire extra letters.
tag<-function(p,l)p+plot_annotation(tag_levels=list(l))
tag <- function(p,l) wrap_elements(full=p)+labs(tag=l)+theme(plot.tag=element_text(face='bold',size=15),plot.tag.position=c(0,1))
main <- wrap_plots(A=tag(pc,'a'),B=tag(pb,'b'),C=tag(pd,'c'),D=tag(pe,'d'),E=tag(pf,'e'),design='AB\nCD\nCE',heights=c(1,1,1))
export(main,'Figure4_assembled_v2',350,360)
export(clean(pStage),'FigS_stage_scores_v2',160,100);export(clean(pFull)+scale_x_discrete(labels=function(z)sub('SAM-','Shoot apex-',z)),'FigS_full_module_correlations_v2',240,115)
# Full dynamic-gene heatmap from existing significant gene set, peak-ordered.
ids<-intersect(dyn$gene[!is.na(dyn$q_value) & dyn$q_value<.01 & dyn$morans_I>0],rownames(s[['SCT']]))
aa<-GetAssayData(s,assay='SCT',layer='data')[ids,spots,drop=FALSE]
bb<-sapply(1:20,function(b)Matrix::rowMeans(aa[,bins==b,drop=FALSE]));zz<-t(scale(t(bb)));zz[!is.finite(zz)]<-0
ord<-order(max.col(zz,ties.method='first'));zz<-zz[ord,,drop=FALSE]
hh<-data.frame(gene=rep(rownames(zz),20),row=rep(seq_len(nrow(zz)),20),bin=rep(1:20,each=nrow(zz)),z=as.vector(zz))
ph<-ggplot(hh,aes(bin,-row,fill=z))+geom_raster()+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-2,2),oob=scales::squish)+theme_minimal(base_size=8)+labs(x='Developmental pseudotime →',y=paste(length(ids),'trajectory-associated genes'),fill='Gene z-score')+theme(axis.text.y=element_blank(),panel.grid=element_blank())
export(ph,'FigS_full_trajectory_gene_heatmap_v2',150,210);write.csv(hh,file.path(out,'FigS_full_trajectory_gene_heatmap_source.csv'),row.names=FALSE)
# Supporting module activity and biological functions.
pmod<-SpatialFeaturePlot(s,features='ME_SAM-M4',images='sam_A_s1',crop=TRUE,pt.size.factor=1.75)+labs(title='Shoot apex-M4',fill='Module eigengene')+theme(legend.position='bottom',plot.title=element_text(size=9),legend.title=element_text(size=6),legend.text=element_text(size=6))+guides(fill=guide_colorbar(title.position='top',barwidth=grid::unit(22,'mm')))
go<-read.csv(file.path(base,'tables/hdWGCNA/SAM_module_GO_fullmodule_2026.7.23.csv'))|>filter(module=='SAM-M4',p.adjust<.05)|>arrange(p.adjust)|>head(12)
go$ratio<-vapply(strsplit(go$GeneRatio,'/'),function(z)as.numeric(z[1])/as.numeric(z[2]),numeric(1))
pgo<-ggplot(go,aes(ratio,reorder(Description,ratio),size=count,color=-log10(p.adjust)))+geom_point()+scale_color_viridis_c()+theme_classic(base_size=8)+labs(x='Gene ratio',y=NULL,color='−log10(FDR)',size='Gene count')
export(pmod+pgo+plot_layout(widths=c(.6,1.4)),'FigS_M4_spatial_activity_and_GO_v2',220,120);write.csv(go,file.path(out,'FigS_M4_GO_source.csv'),row.names=FALSE)
w<-read.csv(file.path(base,'deconvolution/RCTD/tables/SAM_RCTD_multi_celltype_weights.csv'),check.names=FALSE)
ww<-w|>left_join(data.frame(spot=rownames(s@meta.data),domain=s$celltypes),by='spot')|>filter(!is.na(domain))|>pivot_longer(-c(spot,domain),names_to='reference',values_to='weight')|>group_by(domain,reference)|>summarise(weight=mean(weight),.groups='drop')
pw<-ggplot(ww,aes(domain,reference,fill=weight))+geom_tile(color='white')+scale_fill_gradient(low='white',high='#B2182B')+theme_minimal(base_size=8)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())+labs(x='Spatial domain',y='Predicted scRNA-seq cluster',fill='Mean weight')
export(pw,'FigS_all_reference_cluster_predictions_v2',183,170);write.csv(ww,file.path(out,'FigS_all_reference_cluster_predictions_source.csv'),row.names=FALSE)

writeLines(c('Figure 4 v2: a single-cell marker evidence; b spatial score and predicted clusters 14/39; c subcluster map and existing Monocle3 graph; d predefined poplar markers across pseudotime; e core score; f module correlations.','Supplement: stage scores; all module correlations; full significant trajectory-gene heatmap; M4 spatial activity/GO; all reference-cluster predictions.','Pseudotime heatmaps: average SCT log expression in 20 equal-count ordered bins, scaled per gene across bins; color display clipped at ±2. Core-marker panel does not imply significance of each gene.','Existing trajectories/modules retained. Source tables and gene-wise trajectory test results included.','MYB38/SMR1/MATL16 gene assignments follow Giabardo et al. (2026); evidence differs among loci; not every allele has independent promoter validation.'),file.path(out,'README_Figure4_v2.txt'))

})

local({
suppressPackageStartupMessages({library(Seurat);library(qs);library(ggplot2);library(patchwork);library(dplyr);library(tidyr);library(monocle3)})
base <- file.path(ATLAS_WORK_ROOT)
out <- file.path(base,'summary/Fig4_trichome/version2');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x <- readRDS(atlas_input("saved_obj/figure4_plot_inputs.rds"))
list2env(x, environment())
s <- qread(file.path(base,'saved_obj/sobj_shoot_apex_trichome_scores_companion_2026-09-14.qs'))
export <- function(p,n,w=183,h=100){for(ext in c('pdf','svg','png')){tmp<-file.path(tempdir(),paste0(n,'.',ext));ggsave(tmp,p,width=w,height=h,units='mm',device=switch(ext,pdf=cairo_pdf,svg=svglite::svglite,png=ragg::agg_png),dpi=300,limitsize=FALSE,bg='white');stopifnot(file.info(tmp)$size>1000);file.copy(tmp,file.path(out,paste0(n,'.',ext)),overwrite=TRUE)}}
clean <- function(p) p+labs(title=NULL,subtitle=NULL)+theme(plot.title=element_blank(),plot.subtitle=element_blank(),text=element_text(family='Arial',size=8))
pa <- clean(pA)+scale_x_discrete(labels=function(z){z[z=='14: Epidermis']<-'14: Epidermis (putative trichome initials)';z[z=='39: Epidermis']<-'39: Epidermis (putative developing trichomes)';z})
pB <- SpatialFeaturePlot(s,features=c('trichome_core_score','rctd_trichome_initial','rctd_trichome_developing'),images='sam_A_s1',crop=TRUE,combine=FALSE,pt.size.factor=1.75)
for(i in 1:3)pB[[i]]<-pB[[i]]+labs(title=c('Core trichome score','Predicted scRNA-seq cluster 14\nPutative trichome initials','Predicted scRNA-seq cluster 39\nPutative developing trichomes')[i],fill=c('Score','Predicted weight','Predicted weight')[i])+theme(legend.position='bottom',plot.title=element_text(size=8),legend.title=element_text(size=6),legend.text=element_text(size=6))+guides(fill=guide_colorbar(title.position='top',barwidth=grid::unit(20,'mm'),barheight=grid::unit(2,'mm')))
pb <- wrap_plots(pB,nrow=1)
# Existing inferred trajectory; no graph or pseudotime recomputation.
cds <- qread(file.path(base,'trajectory/SAM/monocle3/SAM_upward_leaf_branch_monocle3_cds.qs'))
u <- upward;rownames(u)<-u$spot
common <- intersect(colnames(cds),rownames(u)); stopifnot(length(common)>100); cds <- cds[,common]
colData(cds)$developmental_subcluster <- as.character(u[colnames(cds),'subcluster'])
pcgraph <- plot_cells(cds,color_cells_by='developmental_subcluster',label_cell_groups=FALSE,label_leaves=FALSE,label_roots=FALSE,label_branch_points=FALSE,show_trajectory_graph=TRUE,cell_size=.55)+scale_color_manual(values=subcluster_colors,name='Subcluster')+labs(colour='Subcluster',title='Developmental trajectory')+theme_classic(base_size=8)
s$developmental_subcluster <- factor(u[colnames(s),'subcluster'],levels=names(subcluster_colors))
pcmap <- SpatialDimPlot(s,group.by='developmental_subcluster',images='sam_A_s1',cols=subcluster_colors,crop=TRUE,pt.size.factor=1.75)+labs(title='Spatial subclusters',fill='Subcluster')+theme(legend.position='none',plot.title=element_text(size=9))
pc <- pcmap+pcgraph+plot_layout(widths=c(.7,1.3))
# Marker panel is predefined, not filtered by trajectory-test significance.
mk <- marker_tbl[marker_tbl$include_in_core,]
spots <- u$spot[order(u$monocle3_pseudotime_upward_leaf)]
bins <- dplyr::ntile(seq_along(spots),20)
a <- as.matrix(GetAssayData(s,assay='SCT',layer='data')[mk$spatial_feature,spots,drop=FALSE])
bm <- sapply(1:20,function(b)rowMeans(a[,bins==b,drop=FALSE]));z<-t(scale(t(bm)));z[!is.finite(z)]<-0
hm <- as.data.frame(as.table(z));names(hm)<-c('gene','bin','z');hm$bin<-rep(1:20,each=nrow(z));hm$label<-factor(paste0(mk$gene_name[match(hm$gene,mk$spatial_feature)],' | ',mk$gene_id[match(hm$gene,mk$spatial_feature)]),levels=rev(paste0(mk$gene_name,' | ',mk$gene_id)))
# Full trajectory-significant set; marker labels point to actual rows.
full_src <- file.path(out,'FigS_full_trajectory_gene_heatmap_source.csv')
if (TRUE) {
 dg <- read.csv(file.path(base,'tables/trajectory/SAM_upward_leaf_branch_dynamic_genes.csv'))
 ids <- intersect(dg$gene[!is.na(dg$q_value) & dg$q_value<.01 & dg$morans_I>0],rownames(s[['SCT']]))
 aa <- GetAssayData(s,assay='SCT',layer='data')[ids,spots,drop=FALSE]
 bb <- sapply(1:20,function(b)Matrix::rowMeans(aa[,bins==b,drop=FALSE]))
 zz <- t(scale(t(bb)));zz[!is.finite(zz)]<-0
 zz <- zz[order(max.col(zz,ties.method='first')),,drop=FALSE]
 hh <- data.frame(gene=rep(rownames(zz),20),row=rep(seq_len(nrow(zz)),20),bin=rep(1:20,each=nrow(zz)),z=as.vector(zz))
} else hh <- read.csv(full_src)
ann <- unique(hh[c('gene','row')]) |> inner_join(mk,by=c('gene'='spatial_feature')) |> arrange(row)
stopifnot(nrow(ann)==6)
ann$label <- paste0(ann$gene_name,' | ',ann$gene_id)
ann$label_y <- -ann$row
sep <- max(hh$row)*.065
for(i in seq_len(nrow(ann))[-1]) ann$label_y[i] <- min(ann$label_y[i],ann$label_y[i-1]-sep)
# Shift crowded labels upward if needed, preserving their order and leader lines.
if(min(ann$label_y) < -max(hh$row)) ann$label_y <- ann$label_y + (-max(hh$row)-min(ann$label_y))
pd <- ggplot(hh,aes(bin,-row,fill=z))+geom_raster()+
 geom_segment(data=ann,aes(x=20.6,xend=22,y=-row,yend=label_y),inherit.aes=FALSE,linewidth=.25,color='#555555')+
 geom_text(data=ann,aes(x=22.3,y=label_y,label=label),inherit.aes=FALSE,hjust=0,size=2.6,fontface='italic')+
 scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-2,2),oob=scales::squish)+
 scale_x_continuous(limits=c(.5,41),breaks=c(1,5,10,15,20),expand=c(0,0))+
 scale_y_continuous(expand=c(.005,0))+
 labs(x='Developmental pseudotime → (20 equal-count bins)',y=paste(max(hh$row),'trajectory-significant genes'),fill='Gene z-score')+
 theme_minimal(base_size=8)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.grid=element_blank(),legend.position='bottom')
write.csv(ann[c('gene','gene_name','row','label_y')],file.path(out,'Fig4d_highlighted_gene_rows_v2.csv'),row.names=FALSE)
write.csv(hh,file.path(out,'Fig4d_full_trajectory_heatmap_source_v2.csv'),row.names=FALSE)
write.csv(hm,file.path(out,'Fig4d_marker_pseudotime_source.csv'),row.names=FALSE)
dyn <- read.csv(file.path(base,'tables/trajectory/SAM_upward_leaf_branch_dynamic_genes.csv'))
write.csv(mk|>left_join(dyn,by=c('spatial_feature'='gene')),file.path(out,'Fig4d_marker_trajectory_test_audit.csv'),row.names=FALSE)
pe <- clean(pE)+labs(y='Core trichome score',x='Developmental pseudotime')
pf <- clean(pF)+scale_x_discrete(labels=function(z)sub('SAM-','Shoot apex-',z))+scale_y_discrete(labels=c('Core trichome score'='Core trichome score','Predicted scRNA-seq cluster 14 contribution (trichome initials)'='Predicted scRNA-seq cluster 14\n(putative trichome initials)','Predicted scRNA-seq cluster 39 contribution (developing trichomes)'='Predicted scRNA-seq cluster 39\n(putative developing trichomes)','Meristem-to-primordium pseudotime'='Developmental pseudotime'))


# Match the original supplementary WGCNA module colors, preserving correlations.
module_palette <- read.csv(file.path(base,'tables/hdWGCNA/SAM_hdWGCNA_module_assignments_2026.7.9.csv')) |> distinct(module,color)
module_colors <- setNames(module_palette$color,module_palette$module)
add_module_strip <- function(p) {
 ord <- ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()
 ord <- sub('Shoot apex-','SAM-',ord,fixed=TRUE)
 stopifnot(all(ord %in% names(module_colors)))
 p + annotate('rect',xmin=seq_along(ord)-.5,xmax=seq_along(ord)+.5,
              ymin=.12,ymax=.42,fill=unname(module_colors[ord]),colour='white',linewidth=.25)
}
pf <- add_module_strip(pf)
write.csv(module_palette,file.path(out,'Fig4f_module_color_key_v2.csv'),row.names=FALSE)

ww<-read.csv(file.path(out,'FigS_all_reference_cluster_predictions_source.csv'))
focus<-ww |> filter(reference %in% c('14: Epidermis','39: Epidermis'))
focus$reference<-factor(focus$reference,levels=c('39: Epidermis','14: Epidermis'),labels=c('Cluster 39: putative developing trichomes','Cluster 14: putative trichome initials'))
psummary<-ggplot(focus,aes(domain,reference,fill=weight))+geom_tile(color='white')+scale_fill_gradient(low='white',high='#B2182B')+theme_minimal(base_size=8)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())+labs(x='Spatial domain',y='Predicted scRNA-seq cluster',fill='Mean predicted weight')
pmod<-SpatialFeaturePlot(s,features='ME_SAM-M4',images='sam_A_s1',crop=TRUE,pt.size.factor=1.75)+labs(title='Shoot apex-M4',fill='Module eigengene')+theme(legend.position='bottom',plot.title=element_text(size=9),legend.title=element_text(size=6),legend.text=element_text(size=6))+guides(fill=guide_colorbar(title.position='top',barwidth=grid::unit(22,'mm')))
# Contrast-only display adjustment: retain pink identity; cap at the 95th percentile.
m4_display_limits <- c(min(s$`ME_SAM-M4`,na.rm=TRUE),unname(quantile(s$`ME_SAM-M4`,.95,na.rm=TRUE)))
pmod <- pmod & scale_fill_gradientn(colors=c('grey95',module_colors[['SAM-M4']],'#CE3975'),values=c(0,.45,1),limits=m4_display_limits,oob=scales::squish,name='Module eigengene',na.value='grey90')
write.csv(data.frame(lower=m4_display_limits[1],upper=m4_display_limits[2],upper_rule='95th percentile; higher values share maximum color',high_color='#CE3975'),file.path(out,'Fig4g_display_scale_v2.csv'),row.names=FALSE)
go<-read.csv(file.path(out,'FigS_M4_GO_source.csv'))
pgo<-ggplot(go,aes(ratio,reorder(Description,ratio),size=count,color=-log10(p.adjust)))+geom_point()+scale_color_viridis_c()+theme_classic(base_size=8)+labs(x='Gene ratio',y=NULL,color='−log10(FDR)',size='Gene count')
export(psummary,'Fig4a_predicted_cluster_domain_summary_v2',183,70)
export(pb,'Fig4b_spatial_score_and_predictions_v2',183,85)
export(pc,'Fig4c_subclusters_and_trajectory_v2',183,95)
export(pd,'Fig4d_full_trajectory_heatmap_v2',183,190)
export(pe,'Fig4e_core_score_pseudotime_v2',150,100)
export(pf,'Fig4f_module_correlations_v2',210,85)
export(pmod,'Fig4g_M4_spatial_activity_v2',110,90)
export(pgo,'Fig4h_M4_GO_enrichment_v2',170,120)
tag<-function(p,l)wrap_elements(full=p)+labs(tag=l)+theme(plot.tag=element_text(face='bold',size=15),plot.tag.position=c(0,1))
main<-wrap_plots(A=tag(psummary,'a'),B=tag(pb,'b'),C=tag(pc,'c'),D=tag(pd,'d'),E=tag(pe,'e'),F=tag(pf,'f'),G=tag(pmod,'g'),H=tag(pgo,'h'),design='AABB\nCCEE\nDDFF\nDDGH',heights=c(.8,1,1,1))
export(main,'Figure4_rough_design_v2',360,420)
export(main,'Figure4_assembled_v2',360,420)
s2<-wrap_plots(tag(pa,'a'),tag(clean(pStage),'b'),tag(add_module_strip(clean(pFull)),'c'),ncol=1,heights=c(1,1,1))
export(s2,'Supplementary_Figure_S2_reference_and_scores_v2',240,360)
hubs<-c('PtXaTreH.04G133100.v5.1','PtXaAlbH.04G132100.v5.1')
stopifnot(all(hubs %in% rownames(s[['SCT']])))
sp<-SpatialFeaturePlot(s,features=hubs,images='sam_A_s1',crop=TRUE,combine=FALSE,pt.size.factor=1.75)
for(i in 1:2)sp[[i]]<-sp[[i]]+labs(title=sub('.v5.1','',hubs[i],fixed=TRUE))+theme(plot.title=element_text(size=8),legend.position='bottom')
ex<-as.matrix(GetAssayData(s,assay='SCT',layer='data')[hubs,spots])
hubdf<-data.frame(spot=rep(spots,each=2),gene=rep(hubs,length(spots)),expression=as.vector(ex),pseudotime=rep(u[spots,'monocle3_pseudotime_upward_leaf'],each=2))
hp<-ggplot(hubdf,aes(pseudotime,expression))+geom_point(size=.3,alpha=.3)+geom_smooth(method='loess',span=.7,color='black')+facet_wrap(~gene,ncol=1)+theme_classic(base_size=8)+labs(x='Developmental pseudotime',y='SCT log expression')
export((wrap_plots(sp,nrow=1)/hp),'Supplementary_Figure_S3_M4_hub_genes_v2',183,200)
pw<-ggplot(ww,aes(domain,reference,fill=weight))+geom_tile(color='white')+scale_fill_gradient(low='white',high='#B2182B')+theme_minimal(base_size=8)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())+labs(x='Spatial domain',y='Predicted scRNA-seq cluster',fill='Mean weight')
export(pw,'Supplementary_Figure_S4_all_reference_predictions_v2',183,170)
write.csv(focus,file.path(out,'Fig4a_source_v2.csv'),row.names=FALSE)
write.csv(hubdf,file.path(out,'Supplementary_Figure_S3_source_v2.csv'),row.names=FALSE)

})

suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
base <- file.path(ATLAS_WORK_ROOT)
out <- file.path(base,'summary/Fig4_trichome/version2')
d <- read.csv(file.path(out,'main_figure_module_trait_correlations.csv'))
m <- read.csv(file.path(out,'Fig4f_module_color_key_v2.csv'))
ord <- d |> filter(trait=='Core trichome score') |> arrange(correlation) |> pull(module)
traits <- c('Core trichome score','Predicted scRNA-seq cluster 14 contribution (trichome initials)','Predicted scRNA-seq cluster 39 contribution (developing trichomes)','Meristem-to-primordium pseudotime')
stopifnot(setequal(unique(d$trait),traits),length(ord)==17)
d$x <- match(d$module,ord);d$y <- 5-match(d$trait,traits)
cols <- setNames(m$color,m$module)[ord]
p <- ggplot(d,aes(x,y,fill=correlation))+geom_tile(colour='white',linewidth=.25)+geom_text(aes(label=significance),size=2,fontface='bold')+
 annotate('rect',xmin=seq_along(ord)-.5,xmax=seq_along(ord)+.5,ymin=.05,ymax=.37,fill=unname(cols),colour='white',linewidth=.2)+
 scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-1,1),name='Pearson r')+
 scale_x_continuous(breaks=1:17,labels=sub('SAM-','Shoot apex-',ord),expand=expansion(add=0))+
 scale_y_continuous(breaks=1:4,labels=c('Developmental pseudotime','Predicted scRNA-seq cluster 39\n(putative developing trichomes)','Predicted scRNA-seq cluster 14\n(putative trichome initials)','Core trichome score'),expand=expansion(add=0))+
 coord_fixed(ratio=1,ylim=c(0,4.5),xlim=c(.5,17.5),clip='off')+
 labs(x=NULL,y=NULL)+theme_minimal(base_size=8,base_family='Arial')+
 theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.text.y=element_text(colour='black',lineheight=.95),axis.ticks=element_blank(),legend.position='bottom',legend.title=element_text(size=8),plot.margin=margin(3,4,3,4))+
 guides(fill=guide_colorbar(title.position='left',barwidth=grid::unit(48,'mm'),barheight=grid::unit(2.5,'mm')))
for(ext in c('pdf','svg','png')){
 tmp <- file.path(tempdir(),paste0('Fig4f_module_correlations_compact_v2.',ext))
 ggsave(tmp,p,width=230,height=72,units='mm',device=switch(ext,pdf=cairo_pdf,svg=svglite::svglite,png=ragg::agg_png),dpi=300,bg='white')
 file.copy(tmp,file.path(out,basename(tmp)),overwrite=TRUE)
}

  },
  figure5 = function() {
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
b<-file.path(ATLAS_WORK_ROOT);src<-file.path(b,'summary/Fig5_petiole/two_tissue_development_v2');o<-file.path(src,'Mfuzz_side_comparison_v2');models<-readRDS(file.path(o,'cluster_comparison_models.rds'));dat<-readRDS(file.path(o,'side_resolved_TMM_profiles.rds'));mem<-read.csv(file.path(o,'all_cluster_memberships.csv'));go<-read.csv(file.path(o,'side_DEG_GO_expressed_background_all.csv'));sts<-c('L1','L2','L4','L5','L15');pal<-c(adaxial='#00BFC4',abaxial='#F8766D');ann<-read.csv(file.path(src,'downstream_epidermis_cortex_prioritized_303.csv')) |> select(gene,Functional.Annotation) |> distinct()
export<-function(p,n,w=180,h=110){for(ex in c('pdf','svg','png')){tmp<-file.path(tempdir(),paste0(n,'.',ex));ggsave(tmp,p,width=w,height=h,units='mm',device=switch(ex,pdf=cairo_pdf,svg=svglite::svglite,png=ragg::agg_png),dpi=350,bg='white',limitsize=FALSE);stopifnot(file.info(tmp)$size>100);file.copy(tmp,file.path(o,basename(tmp)),overwrite=TRUE)}}
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold'),legend.position='top',panel.spacing=grid::unit(4,'mm')))
selected<-list();pp<-list();counts<-list();gos<-list();cross<-list()
for(ct in names(dat)){
 k<-if(ct=='Epidermis')3 else 4;md<-models[[paste0(ct,'_Side_DEG_joint')]];f<-md$fits[[as.character(k)]];prefix<-if(ct=='Epidermis')'E' else 'C';mm<-mem |> filter(celltype==ct,analysis=='Side_DEG_joint',.data$k==!!k) |> mutate(panel_cluster=paste0(prefix,sub('G','',cluster)));selected[[ct]]<-mm |> left_join(ann,by='gene');ns<-mm |> group_by(panel_cluster) |> summarise(genes=n(),core_genes=sum(membership>=.5),.groups='drop');counts[[ct]]<-mutate(ns,celltype=ct);labs<-setNames(paste0(ns$panel_cluster,'  (',ns$genes,' genes; ',ns$core_genes,' core)'),ns$panel_cluster)
 z<-as.data.frame(md$x);z$gene<-rownames(md$x);z<-z |> pivot_longer(-gene,names_to='position',values_to='z') |> mutate(side=sub('_.*','',position),stage=factor(sub('^[^_]+_','',position),levels=sts)) |> left_join(mm |> select(gene,panel_cluster,membership),by='gene');c<-as.data.frame(f$centers);c$panel_cluster<-paste0(prefix,1:k);c<-c |> pivot_longer(-panel_cluster,names_to='position',values_to='z') |> mutate(side=sub('_.*','',position),stage=factor(sub('^[^_]+_','',position),levels=sts))
 p<-ggplot(filter(z,membership>=.5),aes(stage,z,color=side,group=interaction(gene,side)))+geom_line(alpha=.16,linewidth=.25)+geom_line(data=c,aes(group=side),linewidth=.95)+geom_point(data=c,aes(group=side),size=1.2)+facet_wrap(~panel_cluster,ncol=if(k==3)3 else 2,labeller=labeller(panel_cluster=labs))+scale_color_manual(values=pal,labels=c(adaxial='Adaxial',abaxial='Abaxial'),name=NULL)+labs(x='Leaf position',y='Standardized expression',title=paste(ct,'adaxial–abaxial DEGs'))
 export(p,paste0('Paired_Mfuzz_',ct,'_k',k,'_v2'),180,if(k==3)75 else 115);pp[[ct]]<-p
 for(cl in ns$panel_cluster){q<-ggplot(filter(z,membership>=.5,panel_cluster==cl),aes(stage,z,color=side,group=interaction(gene,side)))+geom_line(alpha=.18,linewidth=.25)+geom_line(data=filter(c,panel_cluster==cl),aes(group=side),linewidth=1)+scale_color_manual(values=pal,labels=c(adaxial='Adaxial',abaxial='Abaxial'),name=NULL)+labs(x='Leaf position',y='Standardized expression',title=paste(ct,labs[cl]));export(q,paste0('Paired_Mfuzz_',ct,'_',cl,'_v2'),85,70)}
 gg<-go |> filter(celltype==ct,.data$k==!!k) |> mutate(panel_cluster=paste0(prefix,sub('G','',cluster)));gos[[ct]]<-gg

}
sel<-bind_rows(selected);gg<-bind_rows(gos);write.csv(sel,file.path(o,'selected_paired_cluster_genes_annotated.csv'),row.names=FALSE);write.csv(bind_rows(counts),file.path(o,'selected_paired_cluster_sizes.csv'),row.names=FALSE);write.csv(gg,file.path(o,'selected_paired_cluster_GO_all.csv'),row.names=FALSE)
# All significant terms, retaining empty columns for groups without enrichment.
gsig<-gg |> filter(FDR_across_clusters<.05,Count>=3) |> mutate(Description=trimws(Description),panel_cluster=factor(panel_cluster,levels=c('E1','E2','E3','C1','C2','C3','C4')))
q<-ggplot(gsig,aes(panel_cluster,Description))+geom_point(aes(size=Count,color=-log10(FDR_across_clusters)))+scale_x_discrete(drop=FALSE)+scale_color_gradient(low='#8FB9D7',high='#173D70',name='−log10(FDR)')+scale_size_area(max_size=5,name='Gene count')+labs(x='Paired expression group',y=NULL,title='Functions enriched among adaxial–abaxial DEGs')+theme(legend.position='right');export(q,'Paired_Mfuzz_GO_v2',180,90)
# Candidate expression retains amplitudes; points are parent libraries and lines are position means.
can<-data.frame(gene=paste0(c('PtXaAlbH.03G173900','PtXaAlbH.13G121100','PtXaAlbH.07G007900','PtXaTreH.09G107100','PtXaTreH.04G098000','PtXaTreH.12G055400'),'.v5.1'),label=c('Peroxidase\nPtXaAlbH.03G173900','XTH-like\nPtXaAlbH.13G121100','XTH6-like\nPtXaAlbH.07G007900','SAUR-like\nPtXaTreH.09G107100','Expansin A5-like\nPtXaTreH.04G098000','Gibberellin-regulated\nPtXaTreH.12G055400'))
r<-bind_rows(lapply(dat,`[[`,'long')) |> filter(gene%in%can$gene) |> left_join(can,by='gene') |> mutate(stage=factor(stage,levels=sts),side=adaxial_abaxial,label=factor(label,levels=can$label),expression=log2(CPM+1));av<-r |> group_by(gene,celltype,stage,side,label) |> summarise(expression=mean(expression),.groups='drop');write.csv(r,file.path(o,'candidate_side_expression_by_library.csv'),row.names=FALSE)
for(ct in names(dat)){
 q<-ggplot(filter(av,celltype==ct),aes(stage,expression,color=side,group=side))+geom_line(linewidth=.65)+geom_point(data=filter(r,celltype==ct),position=position_dodge(width=.15),size=1.8,alpha=.8)+facet_wrap(~label,ncol=3,scales='free_y',drop=TRUE)+scale_color_manual(values=pal,labels=c(adaxial='Adaxial',abaxial='Abaxial'),name=NULL)+labs(x='Leaf position',y='log2(TMM-normalized CPM + 1)',title=ct);export(q,paste0('Candidate_side_profiles_',ct,'_v2'),180,115)
}
for(g in can$gene){rr<-filter(r,gene==g);aa<-filter(av,gene==g);q<-ggplot(aa,aes(stage,expression,color=side,group=side))+geom_line(linewidth=.65)+geom_point(data=rr,position=position_dodge(width=.15),size=2,alpha=.8)+facet_wrap(~celltype,scales='free_y')+scale_color_manual(values=pal,labels=c(adaxial='Adaxial',abaxial='Abaxial'),name=NULL)+labs(x='Leaf position',y='log2(TMM-normalized CPM + 1)',title=can$label[match(g,can$gene)]);export(q,paste0('Candidate_',sub('.v5.1','',g,fixed=TRUE),'_sides_v2'),130,85)}
# Compact preferred evidence layout, kept as a review design with individual source panels.
q<-(pp$Epidermis / pp$Cortex)+plot_layout(heights=c(1,1.5));export(q,'Paired_Mfuzz_review_layout_v2',185,185)

suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
b<-file.path(ATLAS_WORK_ROOT);o<-file.path(b,'summary/Fig5_petiole/two_tissue_development_v2/Mfuzz_side_comparison_v2');sts<-c('L1','L2','L4','L5','L15');pal<-c(adaxial='#00BFC4',abaxial='#F8766D');export<-function(p,n,w=180,h=110){for(ex in c('pdf','svg','png')){tmp<-file.path(tempdir(),paste0(n,'.',ex));ggsave(tmp,p,width=w,height=h,units='mm',device=switch(ex,pdf=cairo_pdf,svg=svglite::svglite,png=ragg::agg_png),dpi=350,bg='white',limitsize=FALSE);file.copy(tmp,file.path(o,basename(tmp)),overwrite=TRUE)}}
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold'),legend.position='top'))
g<-read.csv(file.path(o,'selected_paired_cluster_GO_all.csv')) |> filter(FDR_across_clusters<.05,Count>=3) |> mutate(panel_cluster=factor(panel_cluster,levels=c('E1','E2','E3','C1','C2','C3','C4')),Description=trimws(Description))
p<-ggplot(g,aes(panel_cluster,Description))+geom_point(aes(size=Count,color=-log10(FDR_across_clusters)))+scale_x_discrete(drop=FALSE)+scale_color_gradient(low='#8FB9D7',high='#173D70',name='−log10(FDR)')+scale_size_area(max_size=5,breaks=c(3,4,5),name='Gene count')+labs(x='Paired expression group',y=NULL,title='Functions enriched among adaxial–abaxial DEGs')+theme(legend.position='right');export(p,'Paired_Mfuzz_GO_v2',180,90)
r<-read.csv(file.path(o,'candidate_side_expression_by_library.csv'));r<-r |> filter((gene=='PtXaAlbH.03G173900.v5.1' & celltype=='Epidermis') | (gene=='PtXaAlbH.07G007900.v5.1' & celltype=='Cortex')) |> mutate(stage=factor(stage,levels=sts),label=ifelse(celltype=='Epidermis','Epidermis: peroxidase\nPtXaAlbH.03G173900','Cortex: XTH6-like\nPtXaAlbH.07G007900'))
av<-r |> group_by(gene,celltype,stage,side,label) |> summarise(expression=mean(expression),.groups='drop');p<-ggplot(av,aes(stage,expression,color=side,group=side))+geom_line(linewidth=.7)+geom_point(data=r,position=position_dodge(width=.13),size=1.8)+facet_wrap(~label,ncol=2,scales='free_y')+scale_color_manual(values=pal,breaks=c('adaxial','abaxial'),labels=c('Adaxial','Abaxial'),name=NULL)+labs(x='Leaf position',y='log2(TMM-normalized CPM + 1)');export(p,'Selected_side_development_candidates_v2',170,80);write.csv(r,file.path(o,'Selected_side_development_candidates_source.csv'),row.names=FALSE)

# Recurrence of the 221 prioritized epidermal/cortical genes; PDF only.
library(ggplot2)
b <- file.path(ATLAS_WORK_ROOT)
p <- file.path(b,'summary/Fig5_petiole/two_tissue_development_v2')
d <- read.csv(file.path(p,'downstream_epidermis_cortex_prioritized_303.csv'))
r <- read.csv(file.path(p,'prioritized_gene_recurrence_two_tissues_v2.csv'))
r$category <- ifelse(r$n_contrasts==1,'One tissue–position combination',ifelse(r$n_celltypes>1 & r$n_positions>1,'Both tissues and multiple positions',ifelse(r$n_celltypes>1,'Both tissues at one position','One tissue at multiple positions')))
r$direction <- vapply(r$gene,function(g){x<-d$side_direction[d$gene==g];n<-table(x);if(length(n)>1 && length(unique(as.numeric(n)))==1)'Mixed (equal counts)' else names(n)[which.max(n)]},character(1))
stopifnot(nrow(r)==221,!anyDuplicated(r$gene))
tab <- as.data.frame(table(r$category,r$direction));names(tab)<-c('Category','Direction','Genes');tab<-tab[tab$Genes>0,]
tab$Category<-factor(tab$Category,levels=rev(c('One tissue–position combination','One tissue at multiple positions','Both tissues at one position','Both tissues and multiple positions')))
g<-ggplot(tab,aes(Genes,Category,fill=Direction))+geom_col(width=.65)+scale_fill_manual(values=c('Adaxial up'='#00BFC4','Abaxial up'='#F8766D','Mixed (equal counts)'='#999999'))+scale_x_continuous(expand=expansion(mult=c(0,.04)))+labs(x='Number of genes',y=NULL,fill='Predominant direction')+theme_classic(base_size=10,base_family='Arial')+theme(legend.position='top',axis.ticks.y=element_blank())
ggsave(file.path(p,'Supplementary_FigS5_recurrence_two_tissues_v2.pdf'),g,width=180,height=83,units='mm',device=cairo_pdf)
write.csv(tab,file.path(p,'Supplementary_FigS5_source_v2.csv'),row.names=FALSE)
if(requireNamespace('pdftools',quietly=TRUE))pdftools::pdf_convert(file.path(p,'Supplementary_FigS5_recurrence_two_tissues_v2.pdf'),format='png',dpi=110,filenames='/tmp/poplar_supp_audit/S5_preview.png',verbose=FALSE)
print(tab)

  }
)
atlas_dispatch(stages)
