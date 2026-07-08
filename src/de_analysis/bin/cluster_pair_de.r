#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

required_args <- c("stage", "seurat")
missing_args <- required_args[vapply(required_args, function(x) is.null(argvs[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

argvs$cluster_col <- if (is.null(argvs$cluster_col)) "FinalIdents" else argvs$cluster_col
argvs$cluster_1 <- if (is.null(argvs$cluster_1)) "59" else as.character(argvs$cluster_1)
argvs$cluster_2 <- if (is.null(argvs$cluster_2)) "82" else as.character(argvs$cluster_2)
argvs$assay <- if (is.null(argvs$assay)) "RNA" else argvs$assay
argvs$slot <- if (is.null(argvs$slot)) "data" else argvs$slot
argvs$test_use <- if (is.null(argvs$test_use)) "wilcox" else argvs$test_use
argvs$min_cells <- if (is.null(argvs$min_cells)) 3L else as.integer(argvs$min_cells)
argvs$min_pct <- if (is.null(argvs$min_pct)) 0 else as.numeric(argvs$min_pct)
argvs$logfc_threshold <- if (is.null(argvs$logfc_threshold)) 0 else as.numeric(argvs$logfc_threshold)

obj <- readRDS(argvs$seurat)
obj <- UpdateSeuratObject(obj)
if (!argvs$assay %in% Assays(obj)) {
  stop(sprintf("Assay %s not found in Seurat object. Available assays: %s",
               argvs$assay, paste(Assays(obj), collapse = ", ")))
}
if (!argvs$cluster_col %in% colnames(obj@meta.data)) {
  stop("Cluster metadata column not found: ", argvs$cluster_col)
}

DefaultAssay(obj) <- argvs$assay
cluster_ids <- as.character(obj@meta.data[[argvs$cluster_col]])

cells_1 <- colnames(obj)[cluster_ids == argvs$cluster_1]
cells_2 <- colnames(obj)[cluster_ids == argvs$cluster_2]

if (length(cells_1) < argvs$min_cells || length(cells_2) < argvs$min_cells) {
  stop(sprintf(
    "Insufficient cells for %s cluster pair %s vs %s: %s=%d, %s=%d, min_cells=%d",
    argvs$stage, argvs$cluster_1, argvs$cluster_2,
    argvs$cluster_1, length(cells_1), argvs$cluster_2, length(cells_2), argvs$min_cells
  ))
}

cells_keep <- c(cells_1, cells_2)
contrast_obj <- subset(obj, cells = cells_keep)
group_1 <- paste0("cluster_", argvs$cluster_1)
group_2 <- paste0("cluster_", argvs$cluster_2)
contrast_obj[["de_group"]] <- ifelse(Cells(contrast_obj) %in% cells_1, group_1, group_2)
Idents(contrast_obj) <- "de_group"

markers <- FindMarkers(
  contrast_obj,
  ident.1 = group_1,
  ident.2 = group_2,
  assay = argvs$assay,
  slot = argvs$slot,
  test.use = argvs$test_use,
  min.pct = argvs$min_pct,
  logfc.threshold = argvs$logfc_threshold
)

markers_dt <- as.data.table(markers, keep.rownames = "gene")
markers_dt[, `:=`(
  stage = argvs$stage,
  contrast = sprintf("cluster_%s_vs_%s", argvs$cluster_1, argvs$cluster_2),
  cluster_1 = argvs$cluster_1,
  cluster_2 = argvs$cluster_2,
  group_1 = group_1,
  group_2 = group_2,
  n_cells_59 = if (argvs$cluster_1 == "59") length(cells_1) else if (argvs$cluster_2 == "59") length(cells_2) else NA_integer_,
  n_cells_82 = if (argvs$cluster_1 == "82") length(cells_1) else if (argvs$cluster_2 == "82") length(cells_2) else NA_integer_,
  n_cells_cluster_1 = length(cells_1),
  n_cells_cluster_2 = length(cells_2),
  assay = argvs$assay,
  slot = argvs$slot,
  test_use = argvs$test_use,
  q_value = p_val_adj,
  abs_avg_log2FC = abs(avg_log2FC)
)]

markers_dt[, q_sort := fifelse(is.na(q_value), Inf, q_value)]
setorder(markers_dt, q_sort, -abs_avg_log2FC)
markers_dt[, q_sort := NULL]
markers_file <- sprintf(
  "%s_cluster_%s_vs_%s_de_markers.csv",
  argvs$stage, argvs$cluster_1, argvs$cluster_2
)
fwrite(markers_dt, markers_file)
