#!/usr/bin/env Rscript
renv::load("/scratch/ycc520/flyem")

# Functional enrichment analysis for OPC neurons
# Tests associations between temporal origins and functional subsystems
# With proper statistical corrections and visualizations

# Load required libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(openxlsx))

# Source utilities
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
} else {
  source("src/utils.r")
}

# Set up p-value tracking file
pval_file <- "functional_enrichment_pvalues_corrected.xlsx"
# Remove existing file to avoid duplicates
if (file.exists(pval_file)) {
  file.remove(pval_file)
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


# No longer needed - Excel will handle headers automatically

# Collect p-values for later correction
pval_collector <- list()

# Record p-value with metadata
record_pval <- function(
  test_type, subsystem, n_total, n_subsys, test_used, pval,
  statistic = NA_real_, direction = NA_character_,
  null_expectation = NA_real_, alternative = "two.sided",
  fdr_family = "Existing_Fisher", notch_status = "All",
  skip_reason = NA_character_
) {
  pval_collector[[length(pval_collector) + 1]] <<- list(
    test_type = test_type,
    subsystem = subsystem,
    notch_status = notch_status,
    n_total = n_total,
    n_subsys = n_subsys,
    test_used = test_used,
    statistic = statistic,
    direction = direction,
    null_expectation = null_expectation,
    alternative = alternative,
    fdr_family = fdr_family,
    skip_reason = skip_reason,
    pval = pval
  )
}

# Apply BH correction independently within each declared hypothesis family.
correct_collected_pvals <- function() {
  if (length(pval_collector) == 0) return(NULL)

  pval_df <- do.call(rbind.data.frame, pval_collector)
  pval_df$p_fdr <- NA_real_

  for (family in unique(pval_df$fdr_family)) {
    idx <- which(pval_df$fdr_family == family & is.finite(pval_df$pval))
    if (length(idx) > 0) {
      pval_df$p_fdr[idx] <- p.adjust(pval_df$pval[idx], method = "fdr")
    }
  }

  return(pval_df)
}

# Write all p-values with corrections to Excel
write_corrected_pvals_excel <- function(filename, pval_df) {
  if (is.null(pval_df) || nrow(pval_df) == 0) return(NULL)

  excel_df <- pval_df
  excel_df <- excel_df[, c(
    "test_type", "subsystem", "notch_status", "n_total", "n_subsys", "test_used",
    "statistic", "direction", "null_expectation", "alternative",
    "fdr_family", "skip_reason", "pval", "p_fdr"
  )]
  colnames(excel_df) <- c(
    "Test_Type", "Functional_Subsystem", "Notch_Status", "N_Total",
    "N_Subsystem", "Test_Used", "Statistic", "Direction",
    "Null_Expectation", "Alternative", "FDR_Family", "Skip_Reason",
    "P_value_raw", "P_value_FDR"
  )
  
  # Split data by test type
  test_types <- unique(excel_df$Test_Type)
  
  # Create workbook
  wb <- createWorkbook()
  
  # Add each test type as a separate sheet
  for (test_type in test_types) {
    sheet_data <- excel_df[excel_df$Test_Type == test_type, ]
    # Remove Test_Type column since it's now the sheet name
    sheet_data <- sheet_data[, !colnames(sheet_data) %in% "Test_Type"]
    
    # Add worksheet
    addWorksheet(wb, sheetName = test_type)
    writeData(wb, sheet = test_type, x = sheet_data, rowNames = FALSE)
  }
  
  # Save workbook
  saveWorkbook(wb, filename, overwrite = TRUE)
  
  return(pval_df)
}

# Cache exact temporal null distributions by window sizes, scores, and
# subsystem membership count.
temporal_null_cache <- new.env(parent = emptyenv())

# Count allocations k_j satisfying sum(k_j) = m and 0 <= k_j <= n_j.
count_feasible_allocations <- function(window_sizes, m) {
  dp <- numeric(m + 1)
  dp[1] <- 1
  for (window_size in window_sizes) {
    next_dp <- numeric(m + 1)
    for (subtotal in 0:m) {
      if (dp[subtotal + 1] == 0) next
      max_add <- min(window_size, m - subtotal)
      for (add in 0:max_add) {
        next_dp[subtotal + add + 1] <-
          next_dp[subtotal + add + 1] + dp[subtotal + 1]
      }
    }
    dp <- next_dp
  }
  return(as.integer(dp[m + 1]))
}

# Expected runs after averaging uniformly over every ordering within each
# tied temporal window. Empty windows are omitted before this function.
tie_averaged_runs <- function(subsystem_counts, window_sizes) {
  proportions <- subsystem_counts / window_sizes
  within_transitions <- sum(
    2 * subsystem_counts * (window_sizes - subsystem_counts) / window_sizes
  )
  between_transitions <- 0
  if (length(window_sizes) > 1) {
    between_transitions <- sum(
      proportions[-length(proportions)] *
        (1 - proportions[-1]) +
        (1 - proportions[-length(proportions)]) *
        proportions[-1]
    )
  }
  return(1 + within_transitions + between_transitions)
}

# Enumerate the exact conditional permutation null. Each feasible vector of
# per-window subsystem counts is weighted by its multivariate-hypergeometric
# probability, conditional on the observed total subsystem size.
get_exact_temporal_null <- function(window_sizes, scores, m) {
  key <- paste(
    paste(window_sizes, collapse = ","),
    paste(scores, collapse = ","),
    m,
    sep = "|"
  )
  if (exists(key, envir = temporal_null_cache, inherits = FALSE)) {
    return(get(key, envir = temporal_null_cache, inherits = FALSE))
  }

  n_total <- sum(window_sizes)
  n_states <- count_feasible_allocations(window_sizes, m)
  probabilities <- numeric(n_states)
  trend_statistics <- numeric(n_states)
  runs_statistics <- numeric(n_states)
  allocation <- integer(length(window_sizes))
  state_index <- 0L

  score_mean <- sum(window_sizes * scores) / n_total
  trend_variance <- (
    m * (n_total - m) / (n_total * (n_total - 1))
  ) * sum(window_sizes * (scores - score_mean)^2)

  enumerate_allocations <- function(window_index, remaining) {
    if (window_index == length(window_sizes)) {
      if (remaining < 0 || remaining > window_sizes[window_index]) return()
      allocation[window_index] <<- remaining
      state_index <<- state_index + 1L

      probabilities[state_index] <<- exp(
        sum(lchoose(window_sizes, allocation)) - lchoose(n_total, m)
      )
      trend_numerator <- sum(
        scores * (allocation - window_sizes * m / n_total)
      )
      trend_statistics[state_index] <<-
        trend_numerator / sqrt(trend_variance)
      runs_statistics[state_index] <<-
        tie_averaged_runs(allocation, window_sizes)
      return()
    }

    remaining_capacity <- sum(window_sizes[(window_index + 1):length(window_sizes)])
    min_count <- max(0, remaining - remaining_capacity)
    max_count <- min(window_sizes[window_index], remaining)
    if (min_count > max_count) return()

    for (count in min_count:max_count) {
      allocation[window_index] <<- count
      enumerate_allocations(window_index + 1, remaining - count)
    }
  }

  enumerate_allocations(1L, m)
  if (state_index != n_states) {
    stop(
      "Exact temporal enumeration generated ", state_index,
      " states; expected ", n_states
    )
  }
  probabilities <- probabilities / sum(probabilities)
  result <- list(
    probabilities = probabilities,
    trend_statistics = trend_statistics,
    runs_statistics = runs_statistics,
    expected_runs = sum(probabilities * runs_statistics)
  )
  assign(key, result, envir = temporal_null_cache)
  return(result)
}

exact_cochran_armitage_test <- function(
    subsystem_counts, window_sizes, scores) {
  m <- sum(subsystem_counts)
  n_total <- sum(window_sizes)
  score_mean <- sum(window_sizes * scores) / n_total
  trend_variance <- (
    m * (n_total - m) / (n_total * (n_total - 1))
  ) * sum(window_sizes * (scores - score_mean)^2)
  observed <- sum(
    scores * (subsystem_counts - window_sizes * m / n_total)
  ) / sqrt(trend_variance)
  null <- get_exact_temporal_null(window_sizes, scores, m)
  tolerance <- sqrt(.Machine$double.eps)
  p_value <- sum(
    null$probabilities[
      abs(null$trend_statistics) >= abs(observed) - tolerance
    ]
  )
  direction <- if (observed < -tolerance) {
    "Early"
  } else if (observed > tolerance) {
    "Late"
  } else {
    "No direction"
  }

  return(list(
    p.value = min(1, p_value),
    statistic = observed,
    direction = direction,
    null.expectation = 0,
    alternative = "two.sided"
  ))
}

exact_tie_averaged_runs_test <- function(subsystem_counts, window_sizes, scores) {
  observed <- tie_averaged_runs(subsystem_counts, window_sizes)
  null <- get_exact_temporal_null(
    window_sizes, scores, sum(subsystem_counts)
  )
  tolerance <- sqrt(.Machine$double.eps)
  p_value <- sum(
    null$probabilities[
      null$runs_statistics <= observed + tolerance
    ]
  )

  return(list(
    p.value = min(1, p_value),
    statistic = observed,
    direction = "Contiguous concentration",
    null.expectation = null$expected_runs,
    alternative = "less"
  ))
}


# Calculate confidence interval for proportion using Wilson score interval
calc_prop_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) return(c(0, 0))
  
  # Use Wilson score interval which is more reliable for small samples
  p <- x/n
  z <- qnorm((1 + conf.level)/2)
  denominator <- 1 + z^2/n
  center <- (p + z^2/(2*n))/denominator
  margin <- z * sqrt(p*(1-p)/n + z^2/(4*n^2))/denominator
  
  ci_lower <- pmax(0, center - margin)
  ci_upper <- pmin(1, center + margin)
  
  return(c(ci_lower, ci_upper))
}

format_pvalue <- function(value) {
  if (length(value) == 0 || !is.finite(value)) return("NA")
  if (value < 0.001) return("<0.001")
  return(sprintf("%.3f", value))
}

format_plot_qvalue <- function(result_row) {
  if (
    is.null(result_row) ||
    !is.na(result_row$skip_reason) ||
    !is.finite(result_row$p_fdr)
  ) {
    return("NT")
  }
  return(format_pvalue(result_row$p_fdr))
}

# ============================================================================
# MODULAR ANALYSIS FUNCTIONS
# ============================================================================

record_skipped_temporal_test <- function(
    test_type, subsys, notch_status, n_total, n_subsys, test_used,
    alternative, fdr_family, skip_reason,
    direction = NA_character_, null_expectation = NA_real_) {
  record_pval(
    test_type, subsys, n_total, n_subsys, test_used, NA_real_,
    statistic = NA_real_,
    direction = direction,
    null_expectation = null_expectation,
    alternative = alternative,
    fdr_family = fdr_family,
    notch_status = notch_status,
    skip_reason = skip_reason
  )
}

run_temporal_tests_by_notch <- function(data, subsys) {
  results <- list()

  for (notch_status in c("Notch Off", "Notch On")) {
    stratum <- data %>% filter(Notch == notch_status)
    n_total <- nrow(stratum)
    n_subsys <- sum(stratum$is_subsystem)
    n_non_subsys <- n_total - n_subsys
    result_key <- if (notch_status == "Notch Off") "notch_off" else "notch_on"
    results[[result_key]] <- list()

    temporal_counts <- stratum %>%
      group_by(temporal_id) %>%
      summarise(
        n_total = n(),
        n_subsys = sum(is_subsystem),
        .groups = "drop"
      ) %>%
      arrange(temporal_id)
    window_sizes <- temporal_counts$n_total
    subsystem_counts <- temporal_counts$n_subsys
    temporal_scores <- temporal_counts$temporal_id

    common_skip_reason <- NA_character_
    if (n_subsys == 0) {
      common_skip_reason <- "No subsystem members in Notch stratum"
    } else if (n_non_subsys == 0) {
      common_skip_reason <- "No non-subsystem members in Notch stratum"
    } else if (length(window_sizes) < 2) {
      common_skip_reason <- "Fewer than two observed temporal windows"
    }

    if (is.na(common_skip_reason)) {
      tbl_temporal <- table(stratum$temporal_id, stratum$is_subsystem)
      fisher_obj <- tryCatch(
        fisher.test(tbl_temporal),
        error = function(e) e
      )
      if (inherits(fisher_obj, "error")) {
        fisher_skip <- paste("Test error:", conditionMessage(fisher_obj))
        record_skipped_temporal_test(
          "Temporal", subsys, notch_status, n_total, n_subsys,
          "fisher", "two.sided", "Temporal", fisher_skip
        )
      } else {
        results[[result_key]]$temporal <- fisher_obj
        record_pval(
          "Temporal", subsys, n_total, n_subsys, "fisher",
          fisher_obj$p.value,
          alternative = "two.sided",
          fdr_family = "Temporal",
          notch_status = notch_status
        )
      }

      trend_obj <- tryCatch(
        exact_cochran_armitage_test(
          subsystem_counts, window_sizes, temporal_scores
        ),
        error = function(e) e
      )
      if (inherits(trend_obj, "error")) {
        trend_skip <- paste("Test error:", conditionMessage(trend_obj))
        record_skipped_temporal_test(
          "Cochran_Armitage", subsys, notch_status, n_total, n_subsys,
          "exact conditional permutation", "two.sided",
          "Cochran_Armitage", trend_skip, null_expectation = 0
        )
      } else {
        results[[result_key]]$cochran_armitage <- trend_obj
        record_pval(
          "Cochran_Armitage", subsys, n_total, n_subsys,
          "exact conditional permutation", trend_obj$p.value,
          statistic = trend_obj$statistic,
          direction = trend_obj$direction,
          null_expectation = trend_obj$null.expectation,
          alternative = trend_obj$alternative,
          fdr_family = "Cochran_Armitage",
          notch_status = notch_status
        )
      }
    } else {
      record_skipped_temporal_test(
        "Temporal", subsys, notch_status, n_total, n_subsys,
        "fisher", "two.sided", "Temporal", common_skip_reason
      )
      record_skipped_temporal_test(
        "Cochran_Armitage", subsys, notch_status, n_total, n_subsys,
        "exact conditional permutation", "two.sided",
        "Cochran_Armitage", common_skip_reason, null_expectation = 0
      )
    }

    runs_skip_reason <- common_skip_reason
    if (is.na(runs_skip_reason) && n_subsys < 2) {
      runs_skip_reason <- "Fewer than two subsystem members"
    }
    if (is.na(runs_skip_reason)) {
      runs_obj <- tryCatch(
        exact_tie_averaged_runs_test(
          subsystem_counts, window_sizes, temporal_scores
        ),
        error = function(e) e
      )
      if (inherits(runs_obj, "error")) {
        runs_skip <- paste("Test error:", conditionMessage(runs_obj))
        record_skipped_temporal_test(
          "Wald_Wolfowitz", subsys, notch_status, n_total, n_subsys,
          "exact tie-averaged permutation", "less",
          "Wald_Wolfowitz", runs_skip,
          direction = "Contiguous concentration"
        )
      } else {
        results[[result_key]]$wald_wolfowitz <- runs_obj
        record_pval(
          "Wald_Wolfowitz", subsys, n_total, n_subsys,
          "exact tie-averaged permutation", runs_obj$p.value,
          statistic = runs_obj$statistic,
          direction = runs_obj$direction,
          null_expectation = runs_obj$null.expectation,
          alternative = runs_obj$alternative,
          fdr_family = "Wald_Wolfowitz",
          notch_status = notch_status
        )
      }
    } else {
      record_skipped_temporal_test(
        "Wald_Wolfowitz", subsys, notch_status, n_total, n_subsys,
        "exact tie-averaged permutation", "less",
        "Wald_Wolfowitz", runs_skip_reason,
        direction = "Contiguous concentration"
      )
    }
  }

  return(results)
}

# Run all statistical tests for a subsystem
run_all_tests <- function(data, subsys, el_cut) {
  results <- list()

  results$temporal_by_notch <- run_temporal_tests_by_notch(data, subsys)
  
  # Broad temporal association test
  data$is_early <- data$temporal_id < el_cut
  tbl_broad <- table(data$is_early, data$is_subsystem)
  if (nrow(tbl_broad) >= 2 && ncol(tbl_broad) >= 2) {
    test_obj <- tryCatch(
      {fisher.test(tbl_broad)}, 
      error = function(e) {list(p.value = NA)}
    )
    results$broad_temporal <- list(
      test = test_obj,
      n_total = sum(tbl_broad),
      n_subsys = sum(tbl_broad[, "TRUE"])
    )
    record_pval(
      "Broad_Temporal", subsys, results$broad_temporal$n_total, 
      results$broad_temporal$n_subsys, "fisher", test_obj$p.value
    )
  }
  
  # Notch association test
  tbl_notch <- table(data$Notch, data$is_subsystem)
  if (nrow(tbl_notch) >= 2 && ncol(tbl_notch) >= 2) {
    test_obj <- tryCatch(
      {fisher.test(tbl_notch)}, 
      error = function(e) {list(p.value = NA)}
    )
    results$notch <- list(
      test = test_obj,
      n_total = sum(tbl_notch),
      n_subsys = sum(tbl_notch[, "TRUE"])
    )
    record_pval(
      "Notch", subsys, results$notch$n_total, 
      results$notch$n_subsys, "fisher", test_obj$p.value
    )
  }
  
  # Broad temporal association test
  tbl_broad_temp <- table(data$broad_temp, data$is_subsystem)
  if (nrow(tbl_broad_temp) >= 2 && ncol(tbl_broad_temp) >= 2) {
    test_obj <- tryCatch(
      {fisher.test(tbl_broad_temp)}, 
      error = function(e) {list(p.value = NA)}
    )
    results$broad_temp <- list(
      test = test_obj,
      n_total = sum(tbl_broad_temp),
      n_subsys = sum(tbl_broad_temp[, "TRUE"])
    )
    record_pval(
      "Broad_Temp", subsys, results$broad_temp$n_total, 
      results$broad_temp$n_subsys, "fisher", test_obj$p.value
    )
  }
  
  return(results)
}

# Prepare data for broad temporal plot
prepare_broad_temp_plot_data <- function(data, subsys) {
  broad_temp_stacked_data <- data %>%
    filter(broad_temp %in% c("Early", "Late")) %>%
    group_by(broad_temp) %>%
    summarise(
      not_in_subsystem = sum(!is_subsystem) / n(),
      in_subsystem = sum(is_subsystem) / n(),
      .groups = 'drop'
    ) %>%
    tidyr::pivot_longer(
      cols = c(not_in_subsystem, in_subsystem),
      names_to = "category",
      values_to = "prop"
    ) %>%
    mutate(
      category = factor(
        category, levels = c("not_in_subsystem", "in_subsystem")
      ),
      fill_color = factor(
        case_when(
          category == "not_in_subsystem" ~ "Not in subsystem",
          TRUE ~ subsys
        ),
        levels = c("Not in subsystem", subsys)
      )
    )
  
  return(broad_temp_stacked_data)
}

# Prepare data for Notch plot
prepare_notch_plot_data <- function(data, subsys) {
  notch_stacked_data <- data %>%
    filter(Notch %in% c("Notch Off", "Notch On")) %>%
    group_by(Notch) %>%
    summarise(
      not_in_subsystem = sum(!is_subsystem) / n(),
      in_subsystem = sum(is_subsystem) / n(),
      .groups = 'drop'
    ) %>%
    tidyr::pivot_longer(
      cols = c(not_in_subsystem, in_subsystem),
      names_to = "category",
      values_to = "prop"
    ) %>%
    mutate(
      category = factor(
        category, levels = c("not_in_subsystem", "in_subsystem")
      ),
      fill_color = factor(
        case_when(
          category == "not_in_subsystem" ~ "Not in subsystem",
          TRUE ~ subsys
        ),
        levels = c("Not in subsystem", subsys)
      )
    )
  
  return(notch_stacked_data)
}

# Prepare Notch-faceted temporal plot data. Missing windows within a Notch
# stratum remain absent so they render as gaps rather than zero prevalence.
prepare_temporal_plot_data <- function(data, subsys) {
  temporal_levels <- data %>%
    distinct(temporal_id, temporal_label) %>%
    arrange(temporal_id)

  stacked_data <- data %>%
    filter(Notch %in% c("Notch Off", "Notch On")) %>%
    group_by(temporal_id, temporal_label, Notch) %>%
    summarise(
      not_in_subsystem = sum(!is_subsystem) / n(),
      in_subsystem = sum(is_subsystem) / n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(not_in_subsystem, in_subsystem),
      names_to = "category",
      values_to = "prop"
    ) %>%
    mutate(
      category = factor(
        category,
        levels = c("not_in_subsystem", "in_subsystem")
      ),
      fill_color = factor(
        case_when(
          grepl("not_in_subsystem", category) ~ "Not in subsystem",
          TRUE ~ subsys
        ),
        levels = c("Not in subsystem", subsys)
      ),
      Notch = factor(Notch, levels = c("Notch Off", "Notch On")),
      temporal_label = factor(
        temporal_label,
        levels = temporal_levels$temporal_label
      )
    )
  
  return(stacked_data)
}

# Prepare combined data for temporal+Notch analysis
prepare_combined_data <- function(data, subsys) {
  prop_data <- data %>%
    filter(Notch %in% c("Notch Off", "Notch On")) %>%
    group_by(temporal_label, Notch) %>%
    summarise(
      n_yes = sum(is_subsystem),
      n_total = n(),
      prop = n_yes / n_total,
      ci_lower = calc_prop_ci(n_yes, n_total)[1],
      ci_upper = calc_prop_ci(n_yes, n_total)[2],
      .groups = 'drop'
    ) %>%
    mutate(
      subsystem = subsys,
      temporal_label = factor(
        temporal_label,
        levels = unique(data$temporal_label)
      ),
      temporal_notch = paste(temporal_label, Notch, sep = "_")
    )
  
  return(prop_data)
}

# Create broad temporal association plot
create_broad_temp_plot <- function(plot_data, test_results, subsys, pval_df = NULL) {
  p_value <- ifelse(
    is.null(test_results$broad_temp$test$p.value), 
    NA, 
    test_results$broad_temp$test$p.value
  )
  
  # Get FDR-corrected p-value if available
  p_fdr <- p_value
  if (!is.null(pval_df)) {
    # Find the FDR-corrected p-value for this specific test
    idx <- which(pval_df$test_type == "Broad_Temp" & pval_df$subsystem == subsys)
    if (length(idx) > 0) {
      p_fdr <- pval_df$p_fdr[idx[1]]
    }
  }
  
  p <- ggplot(plot_data, aes(x = broad_temp, y = prop)) +
    geom_col(
      aes(fill = fill_color),
      position = position_stack(),
      alpha = 0.8
    ) +
    labs(
      title = subsys,
      subtitle = paste0("Fisher's Exact Test (FDR) p = ", 
                       round(p_fdr, 3)),
      x = "Broad Temporal Window",
      y = "Proportion of neurons"
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_fill_manual(values = setNames(
      c("grey90", scale_fill_subsystem()$palette(0)[[subsys]]),
      c("Not in subsystem", subsys)
    )) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(size = 10, color = "grey50"),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    )
  
  return(p)
}

# Create Notch association plot
create_notch_plot <- function(plot_data, test_results, subsys, pval_df = NULL) {
  p_value <- ifelse(
    is.null(test_results$notch$test$p.value), 
    NA, 
    test_results$notch$test$p.value
  )
  
  # Get FDR-corrected p-value if available
  p_fdr <- p_value
  if (!is.null(pval_df)) {
    # Find the FDR-corrected p-value for this specific test
    idx <- which(pval_df$test_type == "Notch" & pval_df$subsystem == subsys)
    if (length(idx) > 0) {
      p_fdr <- pval_df$p_fdr[idx[1]]
    }
  }
  
  p <- ggplot(plot_data, aes(x = Notch, y = prop)) +
    geom_col(
      aes(fill = fill_color),
      position = position_stack(),
      alpha = 0.8
    ) +
    labs(
      title = subsys,
      subtitle = paste0("Fisher's Exact Test (FDR) p = ", 
                       round(p_fdr, 3)),
      x = "Notch Status",
      y = "Proportion of neurons"
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_fill_manual(values = setNames(
      c("grey90", scale_fill_subsystem()$palette(0)[[subsys]]),
      c("Not in subsystem", subsys)
    )) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(size = 10, color = "grey50"),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    )
  
  return(p)
}

make_temporal_plot_subtitle <- function(pval_df, subsys) {
  notch_labels <- c("Notch Off" = "N-", "Notch On" = "N+")
  annotation_lines <- Map(
    function(notch_status, notch_label) {

      get_result <- function(test_type) {
        idx <- which(
          pval_df$test_type == test_type &
            pval_df$subsystem == subsys &
            pval_df$notch_status == notch_status
        )
        if (length(idx) == 0) return(NULL)
        pval_df[idx[1], , drop = FALSE]
      }

      fisher_row <- get_result("Temporal")
      trend_row <- get_result("Cochran_Armitage")
      runs_row <- get_result("Wald_Wolfowitz")
      trend_direction <- if (
        is.null(trend_row) ||
        !is.na(trend_row$skip_reason) ||
        is.na(trend_row$direction)
      ) {
        "NT"
      } else {
        tolower(trend_row$direction)
      }

      paste0(
        notch_label, ": Fisher q=", format_plot_qvalue(fisher_row),
        "; trend q=", format_plot_qvalue(trend_row),
        " (", trend_direction, ")",
        "; clustering q=", format_plot_qvalue(runs_row)
      )
    },
    names(notch_labels),
    unname(notch_labels)
  )

  return(paste(unlist(annotation_lines), collapse = "\n"))
}

# Create temporal plot with separate Notch Off and Notch On facets.
create_temporal_faceted_plot <- function(
    plot_data, subsys, pval_df = NULL) {
  subtitle <- make_temporal_plot_subtitle(pval_df, subsys)

  p <- ggplot(plot_data, aes(x = temporal_label, y = prop)) +
    geom_col(
      aes(fill = fill_color),
      position = position_stack()
    ) +
    facet_wrap(
      ~Notch,
      ncol = 2,
      drop = FALSE,
      labeller = as_labeller(c("Notch Off" = "N-", "Notch On" = "N+"))
    ) +
    labs(
      title = subsys,
      subtitle = subtitle,
      x = "Temporal Origin",
      y = "Proportion of neurons",
      fill = paste("Involved in\n", subsys)
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(
      values = setNames(
        c("grey90", scale_fill_subsystem()$palette(0)[[subsys]]),
        c("Not in subsystem", subsys)
      ),
      labels = c("No", "Yes")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      legend.position = "right",
      axis.title.x = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      strip.text = element_text(size = 13),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  return(p)
}

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

argvs <- commandArgs(trailingOnly = TRUE, asValue = TRUE)

if (interactive()) {
  argvs$anno <- "data/visual_neurons_anno.csv"
}

# Validate input data exists
if (!file.exists(argvs$anno)) {
  stop(
    "Input data file not found: ",
    argvs$anno
  )
}

# Load OPC neuron annotations
opc_flywire <- fread(argvs$anno)

# Validate required columns
required_cols <- c("temporal_id", "temporal_label", "func", "Notch", "broad_temp")
missing_cols <- setdiff(required_cols, colnames(opc_flywire))
if (length(missing_cols) > 0) {
  stop(
    sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# Filter for neurons with temporal information
opc_known_tid <- opc_flywire[
  temporal_id != 0 & Confident_annotation == "Y" & func != "Unannotated" & 
  broad_temp %in% c("Early", "Late")
]

message(
  sprintf("Loaded %d neurons with temporal information", nrow(opc_known_tid))
)

# Define early/late cutoff
el_cut <- 6  # After Ey/Hbn - represents mid-neurogenesis transition

# Number of temporal windows
ntempw <- length(unique(opc_known_tid$temporal_label))

# Excel file will be created when writing results

# Extract functional subsystems
subsystems <- unique(opc_known_tid$func)  

# ============================================================================
# MAIN ANALYSIS LOOP - PROCESS EACH SUBSYSTEM ONCE
# ============================================================================

# Storage for plots and combined data
notch_plots <- list()
broad_temp_plots <- list()
temporal_faceted_plots <- list()
combined_data_all <- data.frame()

# Storage for test results and plot data
all_test_results <- list()
all_plot_data <- list()

# First pass: collect all test results and prepare data
for (subsys in subsystems) {
  message(sprintf("\nProcessing subsystem: %s", subsys))
  
  # Create binary indicator once
  opc_known_tid$is_subsystem <- opc_known_tid$func == subsys
  
  # Run all statistical tests
  test_results <- run_all_tests(opc_known_tid, subsys, el_cut)
  
  # Skip if no valid test results
  if (length(test_results) == 0) {
    message(paste("  Skipping", subsys, "- insufficient data"))
    next
  }
  
  # Store test results for later use
  all_test_results[[subsys]] <- test_results
  
  # Prepare and store plot data
  all_plot_data[[subsys]] <- list()
  
  # Prepare data for Notch plot
  if (!is.null(test_results$notch)) {
    notch_plot_data <- prepare_notch_plot_data(opc_known_tid, subsys)
    if (nrow(notch_plot_data) > 0) {
      all_plot_data[[subsys]]$notch <- notch_plot_data
    }
  }
  
  # Prepare data for broad temporal plot
  if (!is.null(test_results$broad_temp)) {
    broad_temp_plot_data <- prepare_broad_temp_plot_data(opc_known_tid, subsys)
    if (nrow(broad_temp_plot_data) > 0) {
      all_plot_data[[subsys]]$broad_temp <- broad_temp_plot_data
    }
  }
  
  # Prepare data for Notch-faceted temporal plot
  temporal_plot_data <- prepare_temporal_plot_data(opc_known_tid, subsys)
  if (nrow(temporal_plot_data) > 0) {
    all_plot_data[[subsys]]$temporal <- temporal_plot_data
  }
  
  # Prepare combined data for separated plots
  combined_data <- prepare_combined_data(opc_known_tid, subsys)
  if (nrow(combined_data) > 0) {
    combined_data_all <- rbind(combined_data_all, combined_data)
  }
}

# Apply FDR correction to all p-values
pval_df <- correct_collected_pvals()

# Second pass: create plots with FDR-corrected p-values
for (subsys in names(all_test_results)) {
  test_results <- all_test_results[[subsys]]
  plot_data <- all_plot_data[[subsys]]
  
  # Create Notch plot
  if (!is.null(plot_data$notch)) {
    notch_plots[[subsys]] <- create_notch_plot(
      plot_data$notch, test_results, subsys, pval_df
    )
  }
  
  # Create broad temporal plot
  if (!is.null(plot_data$broad_temp)) {
    broad_temp_plots[[subsys]] <- create_broad_temp_plot(
      plot_data$broad_temp, test_results, subsys, pval_df
    )
  }
  
  # Create Notch-faceted temporal plot with FDR-corrected p-values
  if (!is.null(plot_data$temporal)) {
    temporal_faceted_plots[[subsys]] <- create_temporal_faceted_plot(
      plot_data$temporal, subsys, pval_df
    )
  }
}

# ============================================================================
# SAVE ALL PLOTS
# ============================================================================

# Save notch association plots
if (length(notch_plots) > 0) {
  outfn <- "functional_association_notch.pdf"
  wrap_plots(notch_plots, ncol = 3)
  ggsave(
    outfn, 
    width = 15, 
    height = 4 * ceiling(length(notch_plots)/3)
  )
  message(sprintf("Saved Notch plots to: %s", outfn))
}

# Save broad temporal association plots
if (length(broad_temp_plots) > 0) {
  outfn <- "functional_association_broad_temporal.pdf"
  wrap_plots(broad_temp_plots, ncol = 3)
  ggsave(
    outfn, 
    width = 15, 
    height = 4 * ceiling(length(broad_temp_plots)/3)
  )
  message(sprintf("Saved broad temporal plots to: %s", outfn))
}

# Save Notch-faceted temporal plots
if (length(temporal_faceted_plots) > 0) {
  outfn <- "functional_subsystems_temporal_by_notch_faceted.pdf"
  wrap_plots(temporal_faceted_plots, ncol = 1)
  ggsave(
    outfn, 
    width = 16,
    height = 5 * length(temporal_faceted_plots)
  )
  message(sprintf("Saved Notch-faceted temporal plots to: %s", outfn))
}

# Create separated plots for comparison
if (nrow(combined_data_all) > 0) {
  combined_plots <- list()
  
  for (notch_status in c("Notch Off", "Notch On")) {
    subset_data <- combined_data_all[combined_data_all$Notch == notch_status, ]
    
    if (nrow(subset_data) == 0) {
      message(paste("No data for", notch_status))
      next
    }
    
    p <- ggplot(subset_data, aes(x = temporal_label, y = prop)) +
      geom_col(aes(fill = subsystem), position = "stack", alpha = 0.8) +
      labs(
        title = paste0("Functional Subsystems vs Temporal Origin (", notch_status, ")"),
        x = "Temporal Origin",
        y = "Proportion in Functional Subsystem",
        fill = "Functional\nSubsystem"
      ) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_subsystem() +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
      )
    
    combined_plots[[notch_status]] <- p
  }
  
  # Save separated plots
  if (length(combined_plots) > 0) {
    outfn <- "functional_association_combined_separated.pdf"
    if (length(combined_plots) == 2) {
      wrap_plots(combined_plots, ncol = 2) +
        plot_layout(guides = 'collect')
    } else {
      wrap_plots(combined_plots, ncol = 1) +
        plot_layout(guides = 'collect')
    }
    ggsave(outfn, width = 16, height = 8)
    message(sprintf("Saved separated plots to: %s", outfn))
  }
}

# ============================================================================
# WRITE CORRECTED P-VALUES
# ============================================================================

pval_summary <- write_corrected_pvals_excel(pval_file, pval_df)

# Print summary of significant results
message("\n=== SUMMARY OF RESULTS (FDR-corrected) ===")
if (!is.null(pval_summary)) {
  sig_results <- pval_summary[
    is.finite(pval_summary$p_fdr) & pval_summary$p_fdr < 0.05,
  ]
  if (nrow(sig_results) > 0) {
    message("Significant associations after FDR correction:")
    for (i in 1:nrow(sig_results)) {
      detail <- ""
      stratum_detail <- if (sig_results$notch_status[i] == "All") {
        ""
      } else {
        paste0(" [", sig_results$notch_status[i], "]")
      }
      if (sig_results$test_type[i] == "Cochran_Armitage") {
        detail <- paste0(" (", tolower(sig_results$direction[i]), " bias)")
      } else if (sig_results$test_type[i] == "Wald_Wolfowitz") {
        detail <- sprintf(
          " (runs %.3f; null expectation %.3f)",
          sig_results$statistic[i],
          sig_results$null_expectation[i]
        )
      }
      message(sprintf(
        "  %s - %s%s: p(FDR) = %.4f%s",
        sig_results$test_type[i],
        sig_results$subsystem[i],
        stratum_detail,
        sig_results$p_fdr[i],
        detail
      ))
    }
  } else {
    message("No significant associations after FDR correction")
  }
}

message("\nFunctional enrichment analysis complete!")
message("Results saved to ./")
message("Corrected p-values saved to: ", pval_file)
