#!/usr/bin/env Rscript
renv::load("/scratch/ycc520/flyem")
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(openxlsx))

#' Combine selector depth result files into Excel workbook
#' @param pattern Pattern to match result files
#' @param summary_file Output filename for statistical summary
#' @param combined_file Output filename for combined CSV
#' @param excel_file Output filename for Excel file with multiple sheets
combine_selector_results <- function(pattern = "_selector_depth\\.csv$",
                                    summary_file = "selector_depth_summary.txt",
                                    combined_file = "combined_selector_depth.csv",
                                    excel_file = "selector_depth_results.xlsx") {

  # Find all result files
  result_files <- list.files(pattern = pattern, recursive = TRUE, full.names = TRUE)
  # The default pattern also matches the prior combined output. Never feed an
  # aggregate back into itself when rerunning this script in an existing result
  # directory.
  aggregate_names <- unique(c(basename(combined_file), "combined_selector_depth.csv"))
  result_files <- result_files[!basename(result_files) %in% aggregate_names]

  if (length(result_files) == 0) {
    writeLines("No results to combine", summary_file)
    cat("No result files found matching pattern:", pattern, "\n")
    return(invisible(NULL))
  }

  # Combine all results. Individual files contain provisional per-batch FDR
  # values; replace them below after the complete testing family is available.
  all_results <- rbindlist(lapply(result_files, fread), fill = TRUE)

  # Standardize tested rows: represent blank skip reasons as NA
  all_results[, skip_reason := trimws(skip_reason)]
  all_results[skip_reason == "", skip_reason := NA_character_]

  # Recompute FDR across all selectors in each biological testing family so
  # results are independent of workflow batch size and match the CAM analysis.
  all_results[, `:=`(p_value_fdr = NA_real_, significant_fdr = NA)]
  fdr_groups <- c("neuropil", "syn_type", "notch_category")
  all_results[
    !is.na(p_value_exceeds_threshold),
    `:=`(
      p_value_fdr = p.adjust(p_value_exceeds_threshold, method = "fdr"),
      significant_fdr = p.adjust(
        p_value_exceeds_threshold,
        method = "fdr"
      ) < 0.05
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
  ordered_cols <- c(
    intersect(preferred_cols, names(all_results)),
    setdiff(names(all_results), preferred_cols)
  )
  all_results <- all_results[, ..ordered_cols]

  # Write combined CSV
  fwrite(all_results, combined_file)
  cat("Combined CSV written to:", combined_file, "\n")

  # Export to Excel with sheets split by neuropil
  # Drop internal QC columns from report workbook only
  excel_results <- copy(all_results)
  excel_drop_cols <- c(
    "conflict_with_references", "overlap_superficial", "overlap_deep",
    "overlap_any", "batch_id"
  )
  excel_results[, (intersect(excel_drop_cols, colnames(excel_results))) := NULL]
  wb <- createWorkbook()

  for (np in sort(unique(excel_results$neuropil))) {
    # Create sheet name (ME_L, LO_L, etc.)
    sheet_name <- np

    addWorksheet(wb, sheet_name)
    sheet_data <- excel_results[neuropil == np & is.na(skip_reason)]
    writeData(wb, sheet_name, sheet_data)
    if (ncol(sheet_data) > 0) {
      setColWidths(wb, sheet_name, cols = seq_len(ncol(sheet_data)), widths = "auto")
    }
  }

  saveWorkbook(wb, excel_file, overwrite = TRUE)
  cat(
    "Excel file written to:", excel_file, "with",
    length(unique(excel_results$neuropil)), "sheets\n"
  )

  # Generate summary statistics
  sink(summary_file)

  cat("Selector Depth Analysis Summary\n")
  cat("================================\n\n")

  cat("Total analyses:", nrow(all_results), "\n")
  cat("Unique genes:", length(unique(all_results$types_of_interest)), "\n")
  cat("Neuropils:", paste(sort(unique(all_results$neuropil)), collapse = ", "), "\n")
  cat("Syn types:", paste(sort(unique(all_results$syn_type)), collapse = ", "), "\n")

  if ("batch_id" %in% colnames(all_results)) {
    n_batches <- length(unique(all_results$batch_id))
    cat("Batches processed:", n_batches, "\n")
  }

  cat("\n")

  format_pct <- function(n, d) {
    if (d == 0) return("NA")
    sprintf("%.1f%%", 100 * n / d)
  }

  # Skip reason analysis
  cat("=== Skip Reason Analysis ===\n")
  skip_table <- table(all_results$skip_reason, useNA = "no")
  if (length(skip_table) > 0) {
    cat("Skipped genes:\n")
    for (reason in names(skip_table)) {
      cat(sprintf("  %s: %d\n", reason, skip_table[reason]))
    }
  }
  tested_count <- sum(is.na(all_results$skip_reason))
  cat(sprintf("Successfully tested: %d\n\n", tested_count))

  # Per-neuropil, per-Notch summary (collapsed across syn_type)
  cat("=== Per-Neuropil Notch Summary ===\n")
  for (np in sort(unique(all_results$neuropil))) {
    cat(sprintf("\n%s\n", np))
    np_subset <- all_results[neuropil == np]

    for (notch_cat in c("Notch On", "Notch Off")) {
      notch_subset <- np_subset[notch_category == notch_cat]
      n_total <- nrow(notch_subset)
      n_tested <- sum(is.na(notch_subset$skip_reason))
      tested <- notch_subset[is.na(skip_reason)]

      if ("p_value_fdr" %in% colnames(tested)) {
        n_sig_fdr <- sum(tested$p_value_fdr < 0.05, na.rm = TRUE)
      } else if ("significant_fdr" %in% colnames(tested)) {
        n_sig_fdr <- sum(tested$significant_fdr, na.rm = TRUE)
      } else {
        n_sig_fdr <- 0
      }

      cat(sprintf("  %s\n", notch_cat))
      cat(sprintf("    Tested successfully: %d / %d (%s)\n",
                  n_tested, n_total, format_pct(n_tested, n_total)))
      cat(sprintf("    FDR < 0.05 among tested: %d / %d (%s)\n",
                  n_sig_fdr, n_tested, format_pct(n_sig_fdr, n_tested)))
    }
  }
  cat("\n")

  # Per-neuropil, per-syn_type breakdown
  for (np in sort(unique(all_results$neuropil))) {
    for (st in sort(unique(all_results$syn_type))) {
      subset <- all_results[neuropil == np & syn_type == st]

      cat(sprintf("\n=== %s %s ===\n", np, st))
      cat("Total analyses:", nrow(subset), "\n")

      # Add Notch category breakdown
      for (notch_cat in c("Notch On", "Notch Off")) {
        notch_subset <- subset[notch_category == notch_cat]
        cat(sprintf("\n  --- %s ---\n", notch_cat))
        cat(sprintf("  Total genes: %d\n", nrow(notch_subset)))

        tested <- notch_subset[is.na(skip_reason)]
        if (nrow(tested) > 0) {
          cat("  Tested genes:", nrow(tested), "\n")

          # Direction breakdown
          if ("direction" %in% colnames(tested)) {
            direction_table <- table(tested$direction)
            cat("  Direction breakdown:\n")
            for (dir in names(direction_table)) {
              cat(sprintf("    %s: %d\n", dir, direction_table[dir]))
            }
          }

          # Significance analysis
          if ("p_value_exceeds_threshold" %in% colnames(tested)) {
            valid_p <- !is.na(tested$p_value_exceeds_threshold)
            n_valid <- sum(valid_p)

            if (n_valid > 0) {
              n_sig_raw <- sum(tested$p_value_exceeds_threshold < 0.05, na.rm = TRUE)
              cat(sprintf("  Significant (raw p < 0.05): %d / %d\n", n_sig_raw, n_valid))
            }
          }

          if ("p_value_fdr" %in% colnames(tested)) {
            n_sig_fdr <- sum(tested$p_value_fdr < 0.05, na.rm = TRUE)
            cat(sprintf("  Significant (FDR < 0.05): %d / %d (%s)\n",
                        n_sig_fdr, nrow(tested), format_pct(n_sig_fdr, nrow(tested))))
          } else if ("significant_fdr" %in% colnames(tested)) {
            n_sig_fdr <- sum(tested$significant_fdr, na.rm = TRUE)
            cat(sprintf("  Significant (FDR < 0.05): %d / %d (%s)\n",
                        n_sig_fdr, nrow(tested), format_pct(n_sig_fdr, nrow(tested))))
          }

          # Reference conflict detection
          if ("conflict_with_references" %in% colnames(tested)) {
            n_conflicts <- sum(tested$conflict_with_references, na.rm = TRUE)
            if (n_conflicts > 0) {
              cat(sprintf("  WARNING: Reference conflicts: %d genes\n", n_conflicts))
              conflict_genes <- tested[conflict_with_references == TRUE]$types_of_interest
              cat("    Conflicting genes:", paste(head(conflict_genes, 10), collapse = ", "))
              if (n_conflicts > 10) cat(" ...")
              cat("\n")
            }
          }
        } else {
          cat("  No tests passed edge case filters\n")
        }
      }
    }
  }

  # Top significant results (across all neuropils/syn_types)
  cat("\n\n=== Top Significant Results (FDR < 0.05) ===\n")
  if ("significant_fdr" %in% colnames(all_results)) {
    sig_results <- all_results[significant_fdr == TRUE & !is.na(significant_fdr)]

    if (nrow(sig_results) > 0) {
      # Sort by FDR p-value
      sig_results <- sig_results[order(p_value_fdr)]

      for (i in 1:min(20, nrow(sig_results))) {
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
  }

  sink()

  cat("Summary written to:", summary_file, "\n")
}

# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Set defaults
pattern <- if (is.null(argvs$pattern)) "_selector_depth\\.csv$" else argvs$pattern
summary_file <- if (is.null(argvs$summary)) "selector_depth_summary.txt" else argvs$summary
combined_file <- if (is.null(argvs$combined)) "combined_selector_depth.csv" else argvs$combined
excel_file <- if (is.null(argvs$excel)) "selector_depth_results.xlsx" else argvs$excel

# Run the combination
combine_selector_results(pattern, summary_file, combined_file, excel_file)
