#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

depth_pattern <- if (is.null(argvs$depth_pattern)) "_type_depth\\.csv$" else argvs$depth_pattern
depth_out <- if (is.null(argvs$depth_out)) "combined_type_depth.csv" else argvs$depth_out
summary_out <- if (is.null(argvs$summary)) "de_analysis_summary.txt" else argvs$summary
condition_dir <- if (is.null(argvs$condition_dir)) "." else argvs$condition_dir

read_many <- function(files) {
  if (length(files) == 0) return(data.table())
  rbindlist(lapply(files, fread), fill = TRUE)
}

depth_files <- list.files(pattern = depth_pattern, recursive = TRUE, full.names = TRUE)

depth <- read_many(depth_files)

if (nrow(depth) > 0) {
  setorder(depth, neuropil, syn_type, cell_type)
  depth[, skip_reason := trimws(as.character(skip_reason))]
  depth[skip_reason == "", skip_reason := NA_character_]
  depth[, p_value_fdr := NA_real_]
  depth[, significant_fdr := NA]
  depth[!is.na(p_value_exceeds_threshold), p_value_fdr := p.adjust(p_value_exceeds_threshold, method = "fdr"),
    by = .(neuropil, syn_type)
  ]
  depth[!is.na(p_value_fdr), significant_fdr := p_value_fdr < 0.05]

  dir.create(condition_dir, recursive = TRUE, showWarnings = FALSE)
  condition_keys <- unique(depth[, .(neuropil, syn_type)])
  for (idx in seq_len(nrow(condition_keys))) {
    key <- condition_keys[idx]
    condition <- depth[neuropil == key$neuropil & syn_type == key$syn_type]
    condition_name <- sprintf("%s_%s", key$neuropil, key$syn_type)
    out_dir <- file.path(condition_dir, condition_name)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(condition, file.path(out_dir, sprintf("%s_type_depth.csv", condition_name)))
  }
}

fwrite(depth, depth_out)

format_pct <- function(n, d) {
  if (is.na(d) || d == 0) return("NA")
  sprintf("%.1f%%", 100 * n / d)
}

sink(summary_out)
cat("Type Depth Summary\n")
cat("==================\n\n")

cat("Depth files:", length(depth_files), "\n\n")

if (nrow(depth) > 0) {
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
} else {
  cat("No depth rows were produced.\n")
}

sink()

cat("Wrote combined type-depth outputs:\n")
cat("  ", depth_out, "\n", sep = "")
cat("  ", summary_out, "\n", sep = "")
cat("  ", condition_dir, "/\n", sep = "")
