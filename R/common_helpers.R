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
  library(clusterProfiler)
  library(enrichplot)
  library(RcppML)
})

init_project_dirs <- function(root = ".") {
  dirs <- c("QC", "clustering", "marker", "saved_obj", "tables", "trajectory", "go", "nmf")
  walk(file.path(root, dirs), ~ dir.create(.x, showWarnings = FALSE, recursive = TRUE))
}

load_slice <- function(data_dir, slice_name, tissue = NULL, stage = NA_character_, section_type = NA_character_) {
  obj <- Load10X_Spatial(data.dir = data_dir, slice = slice_name)
  obj$orig.ident <- slice_name
  obj$slice_id <- slice_name
  obj$parent_library <- slice_name
  obj$tissue <- tissue %||% slice_name
  obj$stage <- stage
  obj$section_type <- section_type
  obj
}

`%||%` <- function(x, y) if (is.null(x)) y else x

add_basic_qc <- function(sobj) {
  sobj <- PercentageFeatureSet(sobj, pattern = "^ChrMt", col.name = "percent_mt")
  sobj <- PercentageFeatureSet(sobj, pattern = "^ChrPt", col.name = "percent_chlp")
  sobj
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

plot_dbscan_check <- function(seu, eps, minPts, title = NULL, pt.size.factor = 1.8) {
  coords <- GetTissueCoordinates(seu)
  sec <- detect_sections_dbscan(seu, eps = eps, minPts = minPts)
  seu$section_auto <- factor(sec)

  p1 <- SpatialDimPlot(seu, group.by = "section_auto", pt.size.factor = pt.size.factor) +
    ggtitle(ifelse(is.null(title), seu$orig.ident[1], title))

  p2 <- ggplot(coords, aes(x = x, y = y, color = factor(sec[rownames(coords)]))) +
    geom_point(size = 1) +
    scale_y_reverse() +
    coord_equal() +
    theme_classic() +
    labs(color = "DBSCAN", title = "Coordinate view")

  p1 | p2
}

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
  if (!keep_noise) valid_clusters <- valid_clusters[valid_clusters != 0]
  if (length(valid_clusters) == 0) stop("No DBSCAN sections detected for ", base_name)

  centroids <- map_dfr(valid_clusters, function(k) {
    idx <- names(sec)[sec == k]
    tibble(cluster = k,
           x_center = median(coords[idx, xcol]),
           y_center = median(coords[idx, ycol]),
           n_spots = length(idx))
  }) %>%
    arrange(y_center, x_center) %>%
    mutate(section_rank = row_number(),
           new_name = paste0(base_name, "_s", section_rank))

  out <- vector("list", nrow(centroids))
  names(out) <- centroids$new_name

  for (i in seq_len(nrow(centroids))) {
    cl <- centroids$cluster[i]
    new_name <- centroids$new_name[i]
    keep_cells <- names(sec)[sec == cl]
    sub_obj <- subset(seu, cells = keep_cells)
    names(sub_obj@images) <- new_name
    sub_obj$orig.ident <- new_name
    sub_obj$slice_id <- new_name
    sub_obj$parent_library <- base_name
    sub_obj$section_rank <- centroids$section_rank[i]
    sub_obj$section_auto <- paste0("section_", centroids$section_rank[i])
    out[[new_name]] <- sub_obj
  }

  attr(out, "centroids") <- centroids
  out
}

split_libraries_from_params <- function(raw_list, split_params) {
  split_out <- list()
  split_tables <- list()

  for (i in seq_len(nrow(split_params))) {
    nm <- split_params$base_name[i]
    if (!nm %in% names(raw_list)) stop("Missing library in raw_list: ", nm)

    tmp <- split_spatial_sections(
      seu = raw_list[[nm]],
      base_name = nm,
      eps = split_params$eps[i],
      minPts = split_params$minPts[i]
    )
    split_out[names(tmp)] <- tmp
    split_tables[[nm]] <- attr(tmp, "centroids")
  }

  list(objects = split_out, centroids = bind_rows(split_tables, .id = "parent_library"))
}

process_spatial_integration <- function(sobj, resolution = 0.4, dims = 1:30, harmony_var = "orig.ident") {
  sobj <- add_basic_qc(sobj)
  sobj <- SCTransform(sobj, assay = "Spatial", verbose = FALSE)
  sobj <- RunPCA(sobj, assay = "SCT", verbose = FALSE)
  sobj <- RunHarmony(
    object = sobj,
    group.by.vars = harmony_var,
    # assay.use = "SCT",
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

make_marker_heatmap <- function(sobj, marker_df, out_pdf) {
  df_avg <- AverageExpression(sobj, assays = "SCT", layer = "data")$SCT %>% as.matrix()
  keep_var <- apply(df_avg, 1, var) > 0
  z_mtx <- t(scale(t(df_avg[keep_var, , drop = FALSE])))
  z_mtx[!is.finite(z_mtx)] <- 0

  plot_data <- marker_df %>%
    filter(GeneID.v5 %in% rownames(z_mtx)) %>%
    mutate(row_label = paste0(celltype, "_", Name))

  ht_matrix <- z_mtx[plot_data$GeneID.v5, , drop = FALSE]
  rownames(ht_matrix) <- plot_data$row_label

  pdf(out_pdf, height = 18, width = 8)
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
    row_names_gp = grid::gpar(fontsize = 6),
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

plot_annotation_qc <- function(sobj, out_pdf, group_col = "predicted_celltype") {
  p1 <- DimPlot(sobj, reduction = "umap", group.by = "orig.ident", alpha = 0.7) + ggtitle("By section")
  p2 <- DimPlot(sobj, reduction = "umap", group.by = group_col, label = TRUE) + ggtitle("By cluster")
  p3 <- SpatialDimPlot(sobj, group.by = "seurat_clusters", crop = FALSE, ncol = 4)
  p4 <- SpatialDimPlot(sobj, group.by = group_col, crop = FALSE, ncol = 4)
  pdf(out_pdf, width = 16, height = 8)
  print(p1 + p2)
  print(p3)
  print(p4)
  dev.off()
}

plot_known_markers <- function(sobj, marker_genes, out_pdf, images = NULL, ncol = 4, pt.size.factor = 1.8) {
  marker_genes <- intersect(marker_genes, rownames(sobj))
  if (length(marker_genes) == 0) return(invisible(NULL))
  pdf(out_pdf, width = 14, height = max(8, ceiling(length(marker_genes) / ncol) * 3))
  print(
    SpatialFeaturePlot(
      object = sobj,
      features = marker_genes,
      images = images,
      alpha = c(0.1, 1),
      ncol = ncol,
      pt.size.factor = pt.size.factor,
      crop = FALSE
    )
  )
  dev.off()
}

rename_clusters_from_table <- function(
    sobj,
    mapping_df,
    cluster_col = "cluster_label",
    label_col = "celltype_final",
    source_ident_col = "seurat_clusters",
    new_meta_col = "celltypes"
) {
  mapping_df <- mapping_df %>%
    dplyr::mutate(
      !!cluster_col := as.character(.data[[cluster_col]]),
      !!label_col := as.character(.data[[label_col]])
    ) %>%
    dplyr::filter(!is.na(.data[[cluster_col]]), !is.na(.data[[label_col]])) %>%
    dplyr::distinct(.data[[cluster_col]], .keep_all = TRUE)
  
  old_ids <- as.character(sobj[[source_ident_col]][, 1])
  map_idx <- match(old_ids, mapping_df[[cluster_col]])
  new_ids <- mapping_df[[label_col]][map_idx]
  new_ids[is.na(new_ids)] <- old_ids[is.na(new_ids)]
  
  sobj[[new_meta_col]] <- new_ids
  Idents(sobj) <- sobj[[new_meta_col]][, 1]
  sobj
}

assign_stage_from_name <- function(x, levels = c("L1", "L2", "L4", "L5", "L15")) {
  out <- case_when(
    str_detect(x, regex("l15|l16", ignore_case = TRUE)) ~ "L15",
    str_detect(x, regex("l5", ignore_case = TRUE)) ~ "L5",
    str_detect(x, regex("l4", ignore_case = TRUE)) ~ "L4",
    str_detect(x, regex("l2", ignore_case = TRUE)) ~ "L2",
    str_detect(x, regex("l1", ignore_case = TRUE)) ~ "L1",
    TRUE ~ NA_character_
  )
  factor(out, levels = levels)
}

run_stage_de <- function(obj, celltype_name, stage1, stage2, min_cells = 30, assay = "SCT") {
  sub <- subset(obj, subset = celltypes == celltype_name & stage %in% c(stage1, stage2))
  if (ncol(sub) == 0) return(NULL)
  tab <- table(sub$stage)
  if (!all(c(stage1, stage2) %in% names(tab)) || any(tab[c(stage1, stage2)] < min_cells)) return(NULL)
  Idents(sub) <- "stage"
  FindMarkers(sub, ident.1 = stage2, ident.2 = stage1, assay = assay, test.use = "wilcox",
              recorrect_umi = FALSE, logfc.threshold = 0.25, min.pct = 0.1) |>
    tibble::rownames_to_column("gene") |>
    mutate(celltype = celltype_name, stage1 = stage1, stage2 = stage2)
}

run_nmf_programs <- function(sobj, k = 8, assay = "SCT") {
  DefaultAssay(sobj) <- assay
  nmf_features <- intersect(VariableFeatures(sobj), rownames(sobj))
  mat <- GetAssayData(sobj, assay = assay, layer = "data")[nmf_features, ] %>% as.matrix()
  fit <- RcppML::nmf(mat, k = k)
  W <- fit$w
  H <- fit$h
  rownames(W) <- rownames(mat)
  colnames(W) <- paste0("NMF_", seq_len(k))
  rownames(H) <- paste0("NMF_", seq_len(k))
  colnames(H) <- colnames(mat)
  for (i in seq_len(k)) sobj[[paste0("NMF_", i)]] <- H[i, colnames(sobj)]
  list(object = sobj, fit = fit, W = W, H = H, nmf_features = nmf_features)
}

prepare_go_tables <- function(go2gene_file = "genome_files/go2gene_poplar.csv", go2term_file = "genome_files/go2term.csv") {
  go2gene <- read.csv(go2gene_file, header = FALSE)
  go2name <- read.csv(go2term_file, header = TRUE)
  names(go2gene) <- c("ID", "Gene")
  names(go2name) <- c("ID", "Description", "Ontology")
  go2gene$Gene <- ifelse(grepl("\\.v5\\.1$", go2gene$Gene), go2gene$Gene, paste0(go2gene$Gene, ".v5.1"))
  go_anno <- merge(go2gene, go2name, by = "ID")
  go_anno %>% filter(Ontology == "biological_process")
}

run_go_enrich <- function(genes, universe, go_anno_bp, minGSSize = 10, maxGSSize = 500, p_cutoff = 0.05, q_cutoff = 0.2) {
  genes <- intersect(unique(genes), universe)
  if (length(genes) < 5) return(NULL)
  res <- enricher(
    gene = genes,
    universe = universe,
    TERM2GENE = go_anno_bp[, c("ID", "Gene")],
    TERM2NAME = go_anno_bp[, c("ID", "Description")],
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff,
    pAdjustMethod = "BH",
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  as.data.frame(res)
}

run_nmf_gsea <- function(factor_name, W_mat, go_anno_bp, p_cutoff = 0.05) {
  gene_list <- W_mat[, factor_name]
  names(gene_list) <- rownames(W_mat)
  gene_list <- sort(gene_list, decreasing = TRUE)
  GSEA(
    geneList = gene_list,
    TERM2GENE = go_anno_bp[, c("ID", "Gene")],
    TERM2NAME = go_anno_bp[, c("ID", "Description")],
    pvalueCutoff = p_cutoff,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
}

run_endodermis_module <- function(sobj, marker_df, extra_genes = NULL, assay = "SCT", prefix = "endodermis") {
  DefaultAssay(sobj) <- assay
  genes <- marker_df %>% filter(celltype == "Endodermis") %>% pull(GeneID.v5) %>% unique()
  genes <- unique(c(genes, extra_genes))
  genes <- intersect(genes, rownames(sobj))
  if (length(genes) < 3) stop("Too few endodermis genes found in object.")

  sobj <- AddModuleScore(sobj, features = list(genes), name = paste0(prefix, "_score_"), assay = assay, search = FALSE)
  score_col <- paste0(prefix, "_score_1")

  cluster_summary <- sobj@meta.data %>%
    mutate(cluster = as.character(Idents(sobj))) %>%
    group_by(cluster) %>%
    summarise(mean_score = mean(.data[[score_col]], na.rm = TRUE), n = n(), .groups = "drop") %>%
    arrange(desc(mean_score))

  list(object = sobj, genes = genes, score_col = score_col, cluster_summary = cluster_summary)
}


# ===================================================================== #
# Reviewed helper additions for the poplar spatial atlas pipelines
# These functions intentionally keep plotting arguments explicit so that
# SpatialDimPlot/SpatialFeaturePlot layout can be adjusted per tissue.
# ===================================================================== #

make_output_dirs <- function(root = ".") {
  dirs <- c("QC", "section_split", "clustering", "marker", "saved_obj", "tables", "trajectory", "go", "nmf", "figures")
  purrr::walk(file.path(root, dirs), ~ dir.create(.x, showWarnings = FALSE, recursive = TRUE))
}

init_project_dirs <- function(root = ".") {
  make_output_dirs(root)
}

safe_assay_names <- function(sobj) names(sobj@assays)

get_assay_mat <- function(sobj, assay = "SCT", layer = "data") {
  tryCatch(GetAssayData(sobj, assay = assay, layer = layer),
           error = function(e) GetAssayData(sobj, assay = assay, slot = layer))
}

merge_obj_list <- function(obj_list, project) {
  obj_list <- obj_list[!vapply(obj_list, is.null, logical(1))]
  if (length(obj_list) == 0) stop("No objects to merge.")
  if (length(obj_list) == 1) return(obj_list[[1]])
  merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list), project = project)
}

process_spatial_integration_v5 <- function(sobj, resolution = 0.5, dims = 1:30, harmony_vars = "orig.ident",
                                           assay = "Spatial", reduction_save = "harmony") {
  sobj <- add_basic_qc(sobj)
  sobj <- SCTransform(sobj, assay = assay, verbose = FALSE)
  sobj <- RunPCA(sobj, assay = "SCT", verbose = FALSE)
  sobj <- RunHarmony(object = sobj, group.by.vars = harmony_vars, reduction = "pca",
                     reduction.save = reduction_save, plot_convergence = FALSE)
  sobj <- RunUMAP(sobj, reduction = reduction_save, dims = dims)
  sobj <- FindNeighbors(sobj, reduction = reduction_save, dims = dims)
  sobj <- FindClusters(sobj, resolution = resolution)
  sobj
}

process_spatial_integration <- function(sobj, resolution = 0.5, dims = 1:30, harmony_var = "orig.ident") {
  process_spatial_integration_v5(sobj, resolution = resolution, dims = dims, harmony_vars = harmony_var)
}

plot_resolution_sweep <- function(sobj, resolutions, prefix, images_use = NULL,
                                  ncol = 4, pt.size.factor = 4, alpha = 0.8,
                                  crop = FALSE, stroke = 0,
                                  width = 14, height = 12,
                                  out_dir = "clustering") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(images_use)) images_use <- names(sobj@images)[seq_len(min(4, length(names(sobj@images))))]
  pdf(file.path(out_dir, paste0(prefix, "_resolution_optimization.pdf")), width = width, height = height)
  for (res in resolutions) {
    sobj_tmp <- FindClusters(sobj, resolution = res)
    print(DimPlot(sobj_tmp, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
            ggtitle(paste0(prefix, " UMAP | resolution = ", res)))
    for (img in images_use) {
      p <- SpatialDimPlot(sobj_tmp, cells.highlight = CellsByIdentities(sobj_tmp), facet.highlight = TRUE,
                          images = img, ncol = ncol, pt.size.factor = pt.size.factor,
                          alpha = alpha, crop = crop, stroke = stroke) +
        plot_annotation(title = paste0(prefix, " | resolution = ", res, " | ", img))
      print(p)
    }
  }
  dev.off()
}

plot_final_clusters <- function(sobj, prefix, images_use = NULL, group_col = "seurat_clusters",
                                ncol = 4, pt.size.factor = 1.8, crop = FALSE,
                                width = 16, height = 10, out_dir = "clustering") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(images_use)) images_use <- names(sobj@images)[seq_len(min(4, length(names(sobj@images))))]
  p1 <- DimPlot_scCustom(sobj, reduction = "umap", group.by = "orig.ident", alpha = 0.6,
                         ggplot_default_colors = TRUE) + ggtitle("By section")
  p2 <- DimPlot(sobj, reduction = "umap", group.by = group_col, label = TRUE) + ggtitle("By cluster")
  p3 <- SpatialDimPlot(sobj, images = images_use, group.by = group_col, crop = crop,
                       ncol = ncol, pt.size.factor = pt.size.factor)
  pdf(file.path(out_dir, paste0(prefix, "_final_clusters.pdf")), width = width, height = height)
  print(p1 + p2)
  print(p3)
  dev.off()
  pdf(file.path(out_dir, paste0(prefix, "_final_clusters_split_by_cluster.pdf")), width = width, height = max(height, 12))
  for (img in images_use) {
    print(SpatialDimPlot(sobj, cells.highlight = CellsByIdentities(sobj), facet.highlight = TRUE,
                         images = img, ncol = ncol, pt.size.factor = pt.size.factor,
                         alpha = 0.8, crop = crop, stroke = 0) + plot_annotation(title = paste0(prefix, " | ", img)))
  }
  dev.off()
}

clean_marker_table <- function(marker_df) {
  marker_df <- marker_df %>% mutate(across(where(is.character), ~ iconv(.x, from = "", to = "UTF-8", sub = "")))
  if (!"Name" %in% colnames(marker_df)) {
    marker_df$Name <- if ("gene_name" %in% colnames(marker_df)) marker_df$gene_name else marker_df$GeneID.v5
  }
  marker_df
}

make_marker_heatmap <- function(sobj, marker_df, out_pdf, assay = "SCT", layer = "data",
                                height = 18, width = 8, out_table = NULL) {
  marker_df <- clean_marker_table(marker_df)
  df_avg <- AverageExpression(sobj, assays = assay, layer = layer)[[assay]] %>% as.matrix()
  keep_var <- apply(df_avg, 1, var) > 0
  z_mtx <- t(scale(t(df_avg[keep_var, , drop = FALSE])))
  z_mtx[!is.finite(z_mtx)] <- 0
  plot_data <- marker_df %>% filter(GeneID.v5 %in% rownames(z_mtx)) %>%
    mutate(celltype = trimws(as.character(celltype)), Name = trimws(as.character(Name)), GeneID.v5 = trimws(as.character(GeneID.v5))) %>%
    distinct(celltype, GeneID.v5, .keep_all = TRUE) %>%
    arrange(celltype, Name, GeneID.v5) %>%
    mutate(row_label = make.unique(paste0(celltype, " | ", Name, " | ", GeneID.v5)))
  if (nrow(plot_data) == 0) return(invisible(NULL))
  if (!is.null(out_table)) write.csv(plot_data, out_table, row.names = FALSE)
  ht_matrix <- z_mtx[plot_data$GeneID.v5, , drop = FALSE]
  rownames(ht_matrix) <- plot_data$row_label
  pdf(out_pdf, height = height, width = width)
  ht <- ComplexHeatmap::Heatmap(ht_matrix, name = "Z-score", row_split = plot_data$celltype,
                                row_title_rot = 0, row_title_gp = grid::gpar(fontsize = 8),
                                row_gap = grid::unit(3, "mm"), cluster_rows = FALSE,
                                cluster_columns = TRUE, show_row_names = TRUE,
                                row_names_gp = grid::gpar(fontsize = 6),
                                column_names_gp = grid::gpar(fontsize = 9),
                                heatmap_legend_param = list(direction = "horizontal"))
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
  dev.off()
}

score_and_assign_celltypes <- function(sobj, marker_df, assay = "SCT", cluster_col = "seurat_clusters") {
  DefaultAssay(sobj) <- assay
  marker_df <- clean_marker_table(marker_df)
  md <- sobj@meta.data
  md <- md[, !grepl("^score_|^tmpScore_", colnames(md)), drop = FALSE]
  sobj@meta.data <- md
  marker_sets <- marker_df %>%
    distinct(celltype, GeneID.v5) %>%
    filter(GeneID.v5 %in% rownames(sobj)) %>%
    group_by(celltype) %>% summarise(genes = list(unique(GeneID.v5)), .groups = "drop") %>%
    filter(lengths(genes) > 0)
  if (nrow(marker_sets) == 0) stop("No marker genes found in object.")
  score_names <- make.unique(paste0("score_", make.names(marker_sets$celltype)))
  sobj <- AddModuleScore(object = sobj, features = marker_sets$genes, name = "tmpScore_",
                         assay = assay, search = FALSE)
  added_cols <- paste0("tmpScore_", seq_len(nrow(marker_sets)))
  names(sobj@meta.data)[match(added_cols, names(sobj@meta.data))] <- score_names
  cluster_scores <- sobj@meta.data %>% mutate(cluster_label = as.character(.data[[cluster_col]])) %>%
    group_by(cluster_label) %>% summarise(across(all_of(score_names), mean, na.rm = TRUE), .groups = "drop")
  score_long <- cluster_scores %>% tidyr::pivot_longer(cols = all_of(score_names), names_to = "score_name", values_to = "mean_score") %>%
    mutate(celltype = marker_sets$celltype[match(score_name, score_names)]) %>%
    group_by(cluster_label) %>% arrange(desc(mean_score), .by_group = TRUE) %>% mutate(rank = row_number()) %>% ungroup()
  top_assign <- score_long %>% filter(rank == 1) %>%
    select(cluster_label, predicted_celltype = celltype, predicted_score = mean_score)
  sobj$predicted_celltype <- top_assign$predicted_celltype[match(as.character(sobj@meta.data[[cluster_col]]), top_assign$cluster_label)]
  list(object = sobj, cluster_scores = cluster_scores, top_assign = top_assign, score_long = score_long)
}

make_manual_template_from_current_clusters <- function(sobj, ann_res = NULL, file,
                                                       cluster_col = "seurat_clusters") {
  clusters <- sort(unique(as.character(sobj@meta.data[[cluster_col]])))
  if (!is.null(ann_res) && !is.null(ann_res$top_assign)) {
    x <- ann_res$top_assign %>% transmute(cluster_label = as.character(cluster_label),
                                          celltype_predicted = predicted_celltype,
                                          predicted_score,
                                          celltype_final = predicted_celltype)
    missing <- setdiff(clusters, x$cluster_label)
    if (length(missing) > 0) x <- bind_rows(x, tibble(cluster_label = missing, celltype_predicted = missing,
                                                     predicted_score = NA_real_, celltype_final = missing))
    x <- x %>% arrange(cluster_label)
  } else {
    x <- tibble(cluster_label = clusters, celltype_predicted = clusters,
                predicted_score = NA_real_, celltype_final = clusters)
  }
  readr::write_csv(x, file)
  x
}

plot_annotation_qc <- function(sobj, out_pdf, group_col = "predicted_celltype",
                               images = NULL, ncol = 4, pt.size.factor = 1.6,
                               width = 16, height = 10, crop = FALSE) {
  if (is.null(images)) images <- names(sobj@images)[seq_len(min(4, length(names(sobj@images))))]
  p1 <- DimPlot(sobj, reduction = "umap", group.by = "orig.ident", alpha = 0.7) + ggtitle("By section")
  p2 <- DimPlot(sobj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("By cluster")
  pdf(out_pdf, width = width, height = height)
  print(p1 + p2)
  print(SpatialDimPlot(sobj, group.by = "seurat_clusters", images = images, crop = crop,
                       ncol = ncol, pt.size.factor = pt.size.factor))
  if (group_col %in% colnames(sobj@meta.data)) {
    print(SpatialDimPlot(sobj, group.by = group_col, images = images, crop = crop,
                         ncol = ncol, pt.size.factor = pt.size.factor))
  }
  dev.off()
}

plot_known_markers_paged <- function(sobj, marker_genes, out_pdf, images = NULL,
                                     genes_per_page = 6, ncol = 3, pt.size.factor = 1.8,
                                     crop = FALSE, width = 14, height = 10,
                                     min.cutoff = "q05", max.cutoff = "q95") {
  marker_genes <- unique(marker_genes)
  marker_genes <- marker_genes[marker_genes %in% rownames(sobj)]
  if (length(marker_genes) == 0) return(invisible(NULL))
  gene_pages <- split(marker_genes, ceiling(seq_along(marker_genes) / genes_per_page))
  pdf(out_pdf, width = width, height = height)
  for (i in seq_along(gene_pages)) {
    print(SpatialFeaturePlot(sobj, features = gene_pages[[i]], images = images, ncol = ncol,
                             crop = crop, alpha = c(0.1, 1), pt.size.factor = pt.size.factor,
                             min.cutoff = min.cutoff, max.cutoff = max.cutoff) +
            plot_annotation(title = paste0("Markers page ", i, ": ", paste(gene_pages[[i]], collapse = ", "))))
  }
  dev.off()
}

plot_known_marker_dotplot <- function(sobj, marker_df, marker_genes, out_pdf,
                                      cluster_col = "seurat_clusters",
                                      label_cols = c("celltype", "Name", "GeneID.v5"),
                                      width = 12, height = 14, dot.scale = 6,
                                      axis_text_size = 7) {
  marker_df <- clean_marker_table(marker_df)
  marker_genes <- intersect(unique(marker_genes), rownames(sobj))
  label_df <- marker_df %>% filter(GeneID.v5 %in% marker_genes) %>%
    group_by(GeneID.v5) %>%
    summarise(celltype = paste(unique(celltype), collapse = "/"),
              Name = paste(unique(Name), collapse = "/"), .groups = "drop") %>%
    mutate(plot_label = make.unique(paste0(celltype, " | ", Name, " | ", GeneID.v5)))
  features <- label_df$GeneID.v5
  label_map <- stats::setNames(label_df$plot_label, label_df$GeneID.v5)
  pdf(out_pdf, width = width, height = height)
  p <- DotPlot_scCustom(sobj, features = features, group.by = cluster_col,
                        flip_axes = TRUE, dot.scale = dot.scale, x_lab_rotate = TRUE) +
    scale_x_discrete(labels = label_map) + theme(axis.text.y = element_text(size = axis_text_size))
  print(p)
  dev.off()
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

rename_clusters_from_table <- rename_clusters_from_table_direct

load_marker_list_for_scoring <- function() {
  if (file.exists("updated_markerlist_poplar_2.8.26.csv")) {
    clean_marker_table(read.csv("updated_markerlist_poplar_2.8.26.csv"))
  } else {
    mkr <- read.csv("marker_list_poplar_9.27.24.csv")
    if (!"GeneID.v5" %in% colnames(mkr)) mkr$GeneID.v5 <- paste0(mkr$geneID, ".v5.1")
    if (!"celltype" %in% colnames(mkr) && "cell_type" %in% colnames(mkr)) mkr$celltype <- mkr$cell_type
    clean_marker_table(mkr)
  }
}

run_full_nmf <- function(sobj, k = 12, prefix = "TISSUE", assay = "SCT", layer = "data") {
  DefaultAssay(sobj) <- assay
  nmf_features <- intersect(VariableFeatures(sobj), rownames(sobj))
  mat <- get_assay_mat(sobj, assay = assay, layer = layer)[nmf_features, ] |> as.matrix()
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

select_nmf_rank <- function(
    sobj,
    ranks = seq(4, 20, 2),
    n_replicates = 5,
    assay = "SCT",
    layer = "data",
    seed = 123,
    prefix = "TISSUE"
) {
  stopifnot(length(ranks) >= 3, n_replicates >= 2)
  ranks <- sort(unique(as.integer(ranks)))
  nmf_features <- intersect(VariableFeatures(sobj), rownames(sobj))
  mat <- as.matrix(get_assay_mat(sobj, assay = assay, layer = layer)[nmf_features, , drop = FALSE])
  if (any(!is.finite(mat))) stop("NMF rank selection matrix contains non-finite values.")
  if (min(mat) < 0) {
    warning("Negative values found in the NMF input; truncating them to zero.")
    mat[mat < 0] <- 0
  }
  matrix_ss <- sum(mat^2)
  if (matrix_ss == 0) stop("NMF rank selection matrix contains only zeros.")

  normalize_w <- function(w) {
    denom <- sqrt(colSums(w^2))
    denom[denom == 0] <- 1
    sweep(w, 2, denom, "/")
  }
  greedy_component_stability <- function(w_ref, w_other) {
    similarity <- crossprod(normalize_w(w_ref), normalize_w(w_other))
    matched <- numeric(ncol(similarity))
    for (i in seq_len(ncol(similarity))) {
      hit <- which(similarity == max(similarity), arr.ind = TRUE)[1, ]
      matched[i] <- similarity[hit[1], hit[2]]
      similarity[hit[1], ] <- -Inf
      similarity[, hit[2]] <- -Inf
    }
    mean(matched)
  }

  replicate_results <- purrr::map_dfr(ranks, function(k) {
    w_replicates <- vector("list", n_replicates)
    out <- purrr::map_dfr(seq_len(n_replicates), function(rep_id) {
      fit <- RcppML::nmf(
        mat,
        k = k,
        seed = seed + 1000L * k + rep_id,
        verbose = FALSE
      )
      w_replicates[[rep_id]] <<- fit$w
      reconstruction_error <- sqrt(sum((mat - fit$w %*% fit$h)^2) / matrix_ss)
      tibble::tibble(rank = k, replicate = rep_id, reconstruction_error)
    })
    reference_id <- which.min(out$reconstruction_error)
    w_ref <- w_replicates[[reference_id]]
    out$component_stability <- vapply(
      w_replicates,
      function(w) greedy_component_stability(w_ref, w),
      numeric(1)
    )
    out
  })

  rank_summary <- replicate_results %>%
    dplyr::group_by(rank) %>%
    dplyr::summarise(
      mean_reconstruction_error = mean(reconstruction_error),
      sd_reconstruction_error = stats::sd(reconstruction_error),
      mean_component_stability = mean(component_stability),
      sd_component_stability = stats::sd(component_stability),
      .groups = "drop"
    ) %>%
    dplyr::arrange(rank)

  # Locate the reconstruction-error elbow as the greatest distance below the
  # straight line connecting the smallest and largest tested ranks. Stability
  # then downweights elbows that are not reproducible across initializations.
  x <- (rank_summary$rank - min(rank_summary$rank)) /
    (max(rank_summary$rank) - min(rank_summary$rank))
  err_range <- diff(range(rank_summary$mean_reconstruction_error))
  if (err_range == 0) err_range <- 1
  y <- (rank_summary$mean_reconstruction_error - min(rank_summary$mean_reconstruction_error)) / err_range
  # For a decreasing, convex reconstruction-error curve, the elbow lies above
  # the straight endpoint line y = 1 - x after this normalization.
  elbow_strength <- pmax(0, y - (1 - x))
  rank_summary$elbow_strength <- elbow_strength
  rank_summary$selection_score <- elbow_strength * rank_summary$mean_component_stability
  if (all(!is.finite(rank_summary$selection_score)) ||
      max(rank_summary$selection_score, na.rm = TRUE) <= 0) {
    stop(
      "No NMF elbow was detected. Inspect the rank diagnostics and expand the ",
      "tested rank range instead of accepting the smallest tested rank."
    )
  }
  selected_rank <- rank_summary$rank[which.max(rank_summary$selection_score)]
  rank_summary$selected <- rank_summary$rank == selected_rank

  dir.create("nmf", showWarnings = FALSE, recursive = TRUE)
  write.csv(
    replicate_results,
    file.path("nmf", paste0(prefix, "_NMF_rank_selection_replicates.csv")),
    row.names = FALSE
  )
  write.csv(
    rank_summary,
    file.path("nmf", paste0(prefix, "_NMF_rank_selection_summary.csv")),
    row.names = FALSE
  )
  rank_plot <- ggplot2::ggplot(rank_summary, ggplot2::aes(rank)) +
    ggplot2::geom_line(ggplot2::aes(y = mean_reconstruction_error, color = "Reconstruction error")) +
    ggplot2::geom_point(ggplot2::aes(y = mean_reconstruction_error, color = "Reconstruction error"), size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = mean_component_stability, color = "Component stability")) +
    ggplot2::geom_point(ggplot2::aes(y = mean_component_stability, color = "Component stability"), size = 2) +
    ggplot2::geom_vline(xintercept = selected_rank, linetype = 2) +
    ggplot2::scale_color_manual(values = c("Reconstruction error" = "#0072B2", "Component stability" = "#D55E00")) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      y = "Metric value",
      color = NULL,
      title = paste0(prefix, " NMF rank selection; selected k = ", selected_rank),
      subtitle = "Dashed line maximizes reconstruction-error elbow strength × replicate stability"
    )
  ggplot2::ggsave(
    file.path("nmf", paste0(prefix, "_NMF_rank_selection_diagnostics.pdf")),
    rank_plot, width = 8, height = 5
  )

  list(
    selected_rank = selected_rank,
    summary = rank_summary,
    replicates = replicate_results
  )
}

load_go_annotation <- function() {
  go2gene <- read.csv("genome_files/go2gene_poplar.csv", header = FALSE)
  go2name <- read.csv("genome_files/go2term.csv", header = TRUE)
  names(go2gene) <- c("ID", "Gene")
  names(go2name) <- c("ID", "Description", "Ontology")
  go2gene$Gene <- ifelse(grepl("\\.v5\\.1$", go2gene$Gene), go2gene$Gene, paste0(go2gene$Gene, ".v5.1"))
  merge(go2gene, go2name, by = "ID")
}

run_go_enrich <- function(genes, universe, go_anno, ontology = "biological_process") {
  go_anno_use <- go_anno %>% filter(Ontology == ontology)
  genes <- intersect(unique(genes), universe)
  if (length(genes) < 5) return(NULL)
  res <- clusterProfiler::enricher(gene = genes, universe = universe,
                                   TERM2GENE = go_anno_use[, c("ID", "Gene")],
                                   TERM2NAME = go_anno_use[, c("ID", "Description")],
                                   pvalueCutoff = 0.05, pAdjustMethod = "BH")
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  as.data.frame(res)
}

plot_go_top <- function(go_df, group_cols, out_pdf, n_top = 8, width = 12, height = 8) {
  if (is.null(go_df) || nrow(go_df) == 0) return(invisible(NULL))
  plot_df <- go_df %>% group_by(across(all_of(group_cols))) %>% arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
    slice_head(n = n_top) %>% ungroup() %>% mutate(group_label = do.call(paste, c(across(all_of(group_cols)), sep = " | ")))
  p <- ggplot(plot_df, aes(x = group_label, y = Description, size = Count, color = p.adjust)) +
    geom_point() + theme_classic() + labs(x = NULL, y = NULL, size = "Gene count", color = "BH-adjusted p") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out_pdf, p, width = width, height = height)
}
