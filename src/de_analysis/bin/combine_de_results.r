#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

markers_pattern <- if (is.null(argvs$markers_pattern)) "_de_markers\\.csv$" else argvs$markers_pattern
membership_pattern <- if (is.null(argvs$membership_pattern)) "_de_membership\\.csv$" else argvs$membership_pattern
depth_pattern <- if (is.null(argvs$depth_pattern)) "_type_depth\\.csv$" else argvs$depth_pattern
markers_out <- if (is.null(argvs$markers_out)) "combined_de_markers.csv" else argvs$markers_out
membership_out <- if (is.null(argvs$membership_out)) "combined_de_membership.csv" else argvs$membership_out
depth_out <- if (is.null(argvs$depth_out)) "combined_type_depth.csv" else argvs$depth_out
summary_out <- if (is.null(argvs$summary)) "de_analysis_summary.txt" else argvs$summary

read_many <- function(files) {
  if (length(files) == 0) return(data.table())
  rbindlist(lapply(files, fread), fill = TRUE)
}

marker_files <- list.files(pattern = markers_pattern, recursive = TRUE, full.names = TRUE)
membership_files <- list.files(pattern = membership_pattern, recursive = TRUE, full.names = TRUE)
depth_files <- list.files(pattern = depth_pattern, recursive = TRUE, full.names = TRUE)

markers <- read_many(marker_files)
membership <- read_many(membership_files)
depth <- read_many(depth_files)

fwrite(markers, markers_out)
fwrite(membership, membership_out)
fwrite(depth, depth_out)

format_pct <- function(n, d) {
  if (is.na(d) || d == 0) return("NA")
  sprintf("%.1f%%", 100 * n / d)
}

sink(summary_out)
cat("DE Analysis Summary\n")
cat("===================\n\n")

cat("Marker files:", length(marker_files), "\n")
cat("Membership files:", length(membership_files), "\n")
cat("Depth files:", length(depth_files), "\n\n")

if (nrow(depth) > 0) {
  depth[, skip_reason := trimws(as.character(skip_reason))]
  depth[skip_reason == "", skip_reason := NA_character_]
  tested <- depth[is.na(skip_reason)]

  cat("Type Depth\n")
  cat("----------\n")
  cat("Rows:", nrow(depth), "\n")
  cat("Tested:", nrow(tested), "\n")
  cat("Skipped:", nrow(depth) - nrow(tested), "\n")
  if (nrow(tested) > 0) {
    cat("FDR significant:", sum(tested$significant_fdr == TRUE, na.rm = TRUE), "\n")
    cat("Direction counts:\n")
    print(depth[, .N, by = .(neuropil, syn_type, direction)][order(neuropil, syn_type, direction)])
  }
  skip_counts <- depth[!is.na(skip_reason), .N, by = skip_reason][order(-N)]
  if (nrow(skip_counts) > 0) {
    cat("\nSkip reasons:\n")
    print(skip_counts)
  }
  cat("\n")
}

if (nrow(membership) > 0) {
  cat("DE Membership\n")
  cat("-------------\n")
  membership_by <- intersect(
    c("stage", "neuropil", "syn_type", "contrast", "status", "skip_reason"),
    names(membership)
  )
  membership_summary <- membership[, .N, by = membership_by]
  membership_order <- intersect(c("stage", "neuropil", "syn_type", "contrast"), names(membership_summary))
  if (length(membership_order) > 0) setorderv(membership_summary, membership_order)
  print(membership_summary)
  cat("\n")
}

if (nrow(markers) > 0) {
  cat("Markers\n")
  cat("-------\n")
  cat("Rows:", nrow(markers), "\n")
  contrast_cols <- intersect(c("stage", "neuropil", "syn_type", "contrast"), names(markers))
  cat("Contrasts tested:", uniqueN(do.call(paste, c(markers[, ..contrast_cols], sep = "|"))), "\n")
  if ("p_val_adj" %in% names(markers)) {
    sig <- markers[p_val_adj < 0.05]
    cat("Genes with adjusted p < 0.05:", nrow(sig), " / ", nrow(markers),
        " (", format_pct(nrow(sig), nrow(markers)), ")\n", sep = "")
  }
  cat("\nTop rows by adjusted p-value:\n")
  if ("p_val_adj" %in% names(markers)) {
    print(head(markers[order(p_val_adj)], 20))
  } else {
    print(head(markers, 20))
  }
} else {
  cat("Markers\n")
  cat("-------\n")
  cat("No marker rows were produced.\n")
}

sink()

cat("Wrote combined outputs:\n")
cat("  ", markers_out, "\n", sep = "")
cat("  ", membership_out, "\n", sep = "")
cat("  ", depth_out, "\n", sep = "")
cat("  ", summary_out, "\n", sep = "")
