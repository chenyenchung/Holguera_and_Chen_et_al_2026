#!/usr/bin/env Rscript
renv::load("/scratch/ycc520/flyem")
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(openxlsx))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

pattern <- if (is.null(argvs$pattern)) "_cam_depth\\.csv$" else argvs$pattern
summary_file <- if (is.null(argvs$summary)) "cam_depth_summary.txt" else argvs$summary
combined_file <- if (is.null(argvs$combined)) "combined_cam_depth.csv" else argvs$combined
excel_file <- if (is.null(argvs$excel)) "cam_depth_results.xlsx" else argvs$excel

result_files <- list.files(pattern = pattern, recursive = TRUE, full.names = TRUE)

if (length(result_files) == 0) {
  writeLines("No results to combine", summary_file)
  cat("No result files found matching pattern:", pattern, "\n")
  quit(status = 0)
}

all_results <- rbindlist(lapply(result_files, fread), fill = TRUE)

all_results[, skip_reason := trimws(skip_reason)]
all_results[skip_reason == "", skip_reason := NA_character_]

# Individual batches have per-batch FDR from the reused selector engine.
# Recompute after collecting all CAM rows so results are independent of batch size.
all_results[, `:=`(p_value_fdr = NA_real_, significant_fdr = NA)]
fdr_groups <- c("neuropil", "syn_type", "notch_category")
all_results[
  !is.na(p_value_exceeds_threshold),
  `:=`(
    p_value_fdr = p.adjust(p_value_exceeds_threshold, method = "fdr"),
    significant_fdr = p.adjust(p_value_exceeds_threshold, method = "fdr") < 0.05
  ),
  by = fdr_groups
]

preferred_cols <- c(
  "types_of_interest", "notch_category", "observed_bias_ratio",
  "p_value_bias_ratio", "p_value_exceeds_threshold", "p_value_fdr",
  "significant_fdr", "direction", "syn_type", "observed_sup_d",
  "observed_deep_d", "observed_distance_diff", "conflict_with_references",
  "overlap_superficial", "overlap_deep", "overlap_any",
  "observed_delta_thres_base", "observed_delta_thres",
  "bootstrap_distance_diff_median", "bootstrap_distance_diff_lower",
  "bootstrap_distance_diff_upper", "bootstrap_delta_thres_base_median",
  "bootstrap_delta_thres_base_lower", "bootstrap_delta_thres_base_upper",
  "bootstrap_bias_ratio_median", "bootstrap_bias_ratio_lower",
  "bootstrap_bias_ratio_upper", "n_neurons_expressing",
  "n_synapses_expressing", "n_types_expressing", "types_expressing",
  "ref_superficial", "ref_deep", "coefficient", "n_bootstrap",
  "conf_level", "neuropil", "batch_id", "sparse_limit", "analysis_date",
  "skip_reason"
)
ordered_cols <- c(intersect(preferred_cols, names(all_results)),
                  setdiff(names(all_results), preferred_cols))
all_results <- all_results[, ..ordered_cols]

fwrite(all_results, combined_file)
cat("Combined CSV written to:", combined_file, "\n")

excel_results <- copy(all_results)
excel_drop_cols <- c(
  "conflict_with_references", "overlap_superficial", "overlap_deep",
  "overlap_any", "batch_id"
)
excel_results[, (intersect(excel_drop_cols, colnames(excel_results))) := NULL]

wb <- createWorkbook()
for (np in sort(unique(excel_results$neuropil))) {
  addWorksheet(wb, np)
  sheet_data <- excel_results[neuropil == np & is.na(skip_reason)]
  writeData(wb, np, sheet_data)
  if (ncol(sheet_data) > 0) {
    setColWidths(wb, np, cols = seq_len(ncol(sheet_data)), widths = "auto")
  }
}
saveWorkbook(wb, excel_file, overwrite = TRUE)
cat("Excel file written to:", excel_file, "\n")

format_pct <- function(n, d) {
  if (d == 0) return("NA")
  sprintf("%.1f%%", 100 * n / d)
}

sink(summary_file)
cat("CAM Depth Analysis Summary\n")
cat("==========================\n\n")
cat("Total analyses:", nrow(all_results), "\n")
cat("Unique CAM rows:", length(unique(all_results$types_of_interest)), "\n")
cat("Neuropils:", paste(sort(unique(all_results$neuropil)), collapse = ", "), "\n")
cat("Syn types:", paste(sort(unique(all_results$syn_type)), collapse = ", "), "\n")
cat("Batches processed:", length(unique(all_results$batch_id)), "\n\n")

cat("=== Skip Reason Analysis ===\n")
skip_table <- table(all_results$skip_reason, useNA = "no")
if (length(skip_table) > 0) {
  for (reason in names(skip_table)) {
    cat(sprintf("  %s: %d\n", reason, skip_table[reason]))
  }
}
tested_count <- sum(is.na(all_results$skip_reason))
cat(sprintf("Successfully tested: %d\n\n", tested_count))

cat("=== Per-Neuropil Notch Summary ===\n")
for (np in sort(unique(all_results$neuropil))) {
  cat(sprintf("\n%s\n", np))
  np_subset <- all_results[neuropil == np]

  for (notch_cat in c("Notch On", "Notch Off")) {
    notch_subset <- np_subset[notch_category == notch_cat]
    n_total <- nrow(notch_subset)
    n_tested <- sum(is.na(notch_subset$skip_reason))
    tested <- notch_subset[is.na(skip_reason)]
    n_sig_fdr <- sum(tested$p_value_fdr < 0.05, na.rm = TRUE)

    cat(sprintf("  %s\n", notch_cat))
    cat(sprintf("    Tested successfully: %d / %d (%s)\n",
                n_tested, n_total, format_pct(n_tested, n_total)))
    cat(sprintf("    FDR < 0.05 among tested: %d / %d (%s)\n",
                n_sig_fdr, n_tested, format_pct(n_sig_fdr, n_tested)))
  }
}
cat("\n")

for (np in sort(unique(all_results$neuropil))) {
  for (st in sort(unique(all_results$syn_type))) {
    subset <- all_results[neuropil == np & syn_type == st]
    cat(sprintf("\n=== %s %s ===\n", np, st))
    cat("Total analyses:", nrow(subset), "\n")

    for (notch_cat in c("Notch On", "Notch Off")) {
      notch_subset <- subset[notch_category == notch_cat]
      tested <- notch_subset[is.na(skip_reason)]
      cat(sprintf("\n  --- %s ---\n", notch_cat))
      cat(sprintf("  Total CAM rows: %d\n", nrow(notch_subset)))
      cat(sprintf("  Tested CAM rows: %d\n", nrow(tested)))

      if (nrow(tested) > 0) {
        direction_table <- table(tested$direction)
        if (length(direction_table) > 0) {
          cat("  Direction breakdown:\n")
          for (dir in names(direction_table)) {
            cat(sprintf("    %s: %d\n", dir, direction_table[dir]))
          }
        }

        n_sig_raw <- sum(tested$p_value_exceeds_threshold < 0.05, na.rm = TRUE)
        n_sig_fdr <- sum(tested$p_value_fdr < 0.05, na.rm = TRUE)
        cat(sprintf("  Significant (raw p < 0.05): %d / %d\n",
                    n_sig_raw, nrow(tested)))
        cat(sprintf("  Significant (FDR < 0.05): %d / %d (%s)\n",
                    n_sig_fdr, nrow(tested), format_pct(n_sig_fdr, nrow(tested))))

        n_conflicts <- sum(tested$conflict_with_references, na.rm = TRUE)
        if (n_conflicts > 0) {
          conflict_genes <- tested[conflict_with_references == TRUE]$types_of_interest
          cat(sprintf("  WARNING: Reference conflicts: %d CAM rows\n", n_conflicts))
          cat("    Conflicting CAM rows:",
              paste(head(conflict_genes, 10), collapse = ", "))
          if (n_conflicts > 10) cat(" ...")
          cat("\n")
        }
      }
    }
  }
}

cat("\n\n=== Top Significant Results (FDR < 0.05) ===\n")
sig_results <- all_results[significant_fdr == TRUE & !is.na(significant_fdr)]
if (nrow(sig_results) > 0) {
  sig_results <- sig_results[order(p_value_fdr)]
  for (i in seq_len(min(20, nrow(sig_results)))) {
    cat(sprintf("%2d. %s (%s) [%s %s]: direction=%s, FDR p=%.4e, bias_ratio=%.3f\n",
                i,
                sig_results$types_of_interest[i],
                sig_results$notch_category[i],
                sig_results$neuropil[i],
                sig_results$syn_type[i],
                sig_results$direction[i],
                sig_results$p_value_fdr[i],
                sig_results$observed_bias_ratio[i]))
  }
} else {
  cat("No significant results found.\n")
}

sink()
cat("Summary written to:", summary_file, "\n")
