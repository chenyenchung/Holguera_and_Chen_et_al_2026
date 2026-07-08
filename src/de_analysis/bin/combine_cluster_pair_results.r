#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

markers_pattern <- if (is.null(argvs$markers_pattern)) "_cluster_[^/]+_vs_[^/]+_de_markers\\.csv$" else argvs$markers_pattern
markers_out <- if (is.null(argvs$markers_out)) "cluster_59_vs_82_de_table.csv" else argvs$markers_out

marker_files <- list.files(pattern = markers_pattern, recursive = TRUE, full.names = TRUE)
if (length(marker_files) == 0) {
  stop("No cluster-pair marker files found with pattern: ", markers_pattern)
}

markers <- rbindlist(lapply(marker_files, fread), fill = TRUE)
required_cols <- c("gene", "avg_log2FC", "q_value", "abs_avg_log2FC", "stage")
missing_cols <- setdiff(required_cols, names(markers))
if (length(missing_cols) > 0) {
  stop("Missing marker columns: ", paste(missing_cols, collapse = ", "))
}

markers[, q_value := as.numeric(q_value)]
markers[, abs_avg_log2FC := as.numeric(abs_avg_log2FC)]
markers[, q_sort := fifelse(is.na(q_value), Inf, q_value)]
setorder(markers, q_sort, -abs_avg_log2FC)
markers[, q_sort := NULL]

fwrite(markers, markers_out)
cat("Wrote cluster-pair DE table:\n")
cat("  ", markers_out, "\n", sep = "")
