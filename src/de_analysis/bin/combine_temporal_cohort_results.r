#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

markers_pattern <- if (is.null(argvs$markers_pattern)) "_temporal_cohort_de_markers\\.csv$" else argvs$markers_pattern
membership_pattern <- if (is.null(argvs$membership_pattern)) "_temporal_cohort_de_membership\\.csv$" else argvs$membership_pattern
cam_file <- if (is.null(argvs$cam)) file.path(repo_root, "data/P15_CAM.csv") else argvs$cam
q_threshold <- if (is.null(argvs$q_threshold)) 0.05 else as.numeric(argvs$q_threshold)
top_n <- if (is.null(argvs$top_n)) 20L else as.integer(argvs$top_n)

markers_out <- if (is.null(argvs$markers_out)) "combined_temporal_cohort_de_markers.csv" else argvs$markers_out
membership_out <- if (is.null(argvs$membership_out)) "combined_temporal_cohort_de_membership.csv" else argvs$membership_out
candidates_out <- if (is.null(argvs$candidates_out)) "temporal_cohort_cam_candidates.csv" else argvs$candidates_out
summary_out <- if (is.null(argvs$summary)) "temporal_cohort_de_summary.txt" else argvs$summary

read_many <- function(files) {
  if (length(files) == 0) return(data.table())
  rbindlist(lapply(files, fread), fill = TRUE)
}

marker_files <- list.files(pattern = markers_pattern, recursive = TRUE, full.names = TRUE)
membership_files <- list.files(pattern = membership_pattern, recursive = TRUE, full.names = TRUE)

markers <- read_many(marker_files)
membership <- read_many(membership_files)
cam_genes <- fread(cam_file, select = 1)[[1]]
cam_genes <- unique(as.character(cam_genes))

if (nrow(markers) > 0) {
  markers[, is_cam := gene %in% cam_genes]
  markers[, candidate_class := fifelse(is_cam, "CAM", "")]
  markers[, enrichment_direction := fifelse(
    avg_log2FC > 0, "Early",
    fifelse(avg_log2FC < 0, "Late", "none")
  )]
}

fwrite(markers, markers_out)
fwrite(membership, membership_out)

if (nrow(markers) > 0) {
  candidates <- markers[
    is_cam == TRUE & p_val_adj < q_threshold & enrichment_direction %in% c("Early", "Late")
  ]
  candidates[, abs_avg_log2FC := abs(avg_log2FC)]
  setorder(candidates, stage, contrast, enrichment_direction, p_val_adj, -abs_avg_log2FC)

  candidates <- candidates[
    ,
    head(.SD, top_n),
    by = .(stage, contrast, enrichment_direction)
  ]
  candidates[, abs_avg_log2FC := NULL]
} else {
  candidates <- data.table()
}
fwrite(candidates, candidates_out)

format_pct <- function(n, d) {
  if (is.na(d) || d == 0) return("NA")
  sprintf("%.1f%%", 100 * n / d)
}

sink(summary_out)
cat("Temporal Cohort DE Summary\n")
cat("==========================\n\n")
cat("Marker files:", length(marker_files), "\n")
cat("Membership files:", length(membership_files), "\n")
cat("CAM genes:", length(cam_genes), "\n\n")

if (nrow(membership) > 0) {
  cat("Membership\n")
  cat("----------\n")
  membership_by <- intersect(
    c("stage", "contrast", "status", "skip_reason"),
    names(membership)
  )
  print(membership[, .N, by = membership_by][order(stage, contrast, status)])
  cat("\n")
}

if (nrow(markers) > 0) {
  cat("Markers\n")
  cat("-------\n")
  cat("Rows:", nrow(markers), "\n")
  contrast_cols <- intersect(c("stage", "contrast"), names(markers))
  cat("Contrasts tested:", uniqueN(do.call(paste, c(markers[, ..contrast_cols], sep = "|"))), "\n")
  sig <- markers[p_val_adj < q_threshold]
  sig_cam <- sig[is_cam == TRUE]
  cat("Genes with adjusted p < ", q_threshold, ": ", nrow(sig), " / ", nrow(markers),
      " (", format_pct(nrow(sig), nrow(markers)), ")\n", sep = "")
  cat("Significant CAM candidates:", nrow(sig_cam), "\n\n")
  cat("Significant CAM candidates by contrast:\n")
  print(sig_cam[, .N, by = .(stage, contrast, enrichment_direction)][
    order(stage, contrast, enrichment_direction)
  ])
  cat("\nTop CAM candidate rows:\n")
  print(head(candidates[order(stage, contrast, p_val_adj)], 40))
} else {
  cat("Markers\n")
  cat("-------\n")
  cat("No marker rows were produced.\n")
}
sink()

cat("Wrote combined temporal cohort outputs:\n")
cat("  ", markers_out, "\n", sep = "")
cat("  ", membership_out, "\n", sep = "")
cat("  ", candidates_out, "\n", sep = "")
cat("  ", summary_out, "\n", sep = "")
