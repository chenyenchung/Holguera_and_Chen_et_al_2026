#!/usr/bin/env Rscript
renv::load("/scratch/ycc520/flyem")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(Rcpp))



#' Compute pairwise quantile Wasserstein distances between groups
#' @param data Data frame with depth and group columns
#' @param sparse_limit Minimum number of observations per group
#' @param n_quantiles Number of quantiles to use (default 100)
#' @param n_bootstrap Number of bootstrap iterations (default 1000)
#' @param conf_int Confidence interval percentage (default 95)
#' @return Data frame with Wasserstein distance results
quantile_wasserstein_dist <- function(data, sparse_limit = 100, 
                                      n_quantiles = 1000, n_bootstrap = 1000,
                                      conf_int = 95) {
  # Filter groups with sufficient data
  group_counts <- table(data$group)
  valid_groups <- names(group_counts)[group_counts > sparse_limit]
  
  if (length(valid_groups) < 2) {
    warning("Less than 2 groups have sufficient data for testing")
    return(data.frame(
      group1 = character(0), group2 = character(0), 
      wasserstein_median = numeric(0), wasserstein_lower = numeric(0),
      wasserstein_upper = numeric(0), conf_level = numeric(0)
    ))
  }
  
  # Use C++ implementation for all comparisons
  group_data <- lapply(valid_groups, function(g) {
    data[data$group == g, "depth"][[1]]
  })
  
  cpp_results <- wasserstein_dist_cpp(group_data, valid_groups,
                                      n_quantiles, n_bootstrap, conf_int)
  
  cpp_results$conf_level <- conf_int
  
  return(cpp_results)
}


#' Perform subsampled Kolmogorov-Smirnov test
#' @param sample1 Numeric vector of first sample
#' @param sample2 Numeric vector of second sample
#' @param subsample_size Number of points to subsample from each group (default 1000)
#' @param n_iterations Number of iterations to perform (default 1000)
#' @return List with raw p-value statistics
subsampled_ks_test <- function(sample1, sample2, subsample_size = 1000, 
                               n_iterations = 1000) {
  
  # Input validation
  if (length(sample1) < 10 || length(sample2) < 10) {
    stop("Samples must have at least 10 observations each")
  }
  
  if (subsample_size > min(length(sample1), length(sample2))) {
    subsample_size <- min(length(sample1), length(sample2))
    warning("Subsample size reduced to minimum group size: ", subsample_size)
  }
  
  # Perform iterations
  raw_pvals <- numeric(n_iterations)
  
  for (i in 1:n_iterations) {
    # Subsample without replacement (or use full sample if smaller than subsample_size)
    if (length(sample1) <= subsample_size) {
      sub1 <- sample1
    } else {
      sub1 <- sample(sample1, size = subsample_size, replace = FALSE)
    }
    
    if (length(sample2) <= subsample_size) {
      sub2 <- sample2
    } else {
      sub2 <- sample(sample2, size = subsample_size, replace = FALSE)
    }
    
    # Perform KS test
    ks_result <- ks.test(sub1, sub2)
    raw_pvals[i] <- ks_result$p.value
  }
  
  # Calculate summary statistics for raw p-values
  raw_stats <- list(
    median = median(raw_pvals, na.rm = TRUE),
    q25 = quantile(raw_pvals, 0.25, na.rm = TRUE),
    q75 = quantile(raw_pvals, 0.75, na.rm = TRUE)
  )
  
  return(list(
    raw_pvals = raw_stats,
    subsample_size = subsample_size,
    n_iterations = n_iterations
  ))
}

#' Perform pairwise subsampled KS tests between groups
#' @param data Data frame with depth and group columns
#' @param sparse_limit Minimum number of observations per group
#' @param subsample_size Number of points to subsample from each group (default 1000)
#' @param n_iterations Number of iterations to perform (default 1000)
#' @param correction_method Correction method to apply to pairwise comparisons
#' @return Data frame with KS test results
perform_subsampled_ks_tests <- function(data, sparse_limit = 100,
                                        subsample_size = 1000, n_iterations = 1000,
                                        correction_method = "fdr") {
  
  # Filter groups with sufficient data
  group_counts <- table(data$group)
  valid_groups <- names(group_counts)[group_counts > sparse_limit]
  
  if (length(valid_groups) < 2) {
    warning("Less than 2 groups have sufficient data for KS testing")
    return(data.frame())
  }
  
  # Generate all pairwise combinations
  combinations <- combn(valid_groups, 2, simplify = FALSE)
  
  # Initialize results list
  results_list <- list()
  
  for (i in seq_along(combinations)) {
    group1 <- combinations[[i]][1]
    group2 <- combinations[[i]][2]
    
    sample1 <- data[data$group == group1, "depth"][[1]]
    sample2 <- data[data$group == group2, "depth"][[1]]
    
    # Perform subsampled KS test
    ks_result <- subsampled_ks_test(sample1, sample2, subsample_size, 
                                    n_iterations)
    
    # Create result row
    result_row <- data.frame(
      group1 = group1,
      group2 = group2,
      ks_pval_raw_median = ks_result$raw_pvals$median,
      ks_pval_raw_q25 = ks_result$raw_pvals$q25,
      ks_pval_raw_q75 = ks_result$raw_pvals$q75,
      ks_subsample_size = ks_result$subsample_size,
      ks_n_iterations = ks_result$n_iterations,
      stringsAsFactors = FALSE
    )
    
    results_list[[paste(group1, group2, sep = "_")]] <- result_row
  }
  
  if (length(results_list) > 0) {
    # Combine all results
    combined_results <- do.call(rbind, results_list)
    
    # Apply multiple testing correction to all summary statistics across pairwise comparisons
    corrected_median <- p.adjust(combined_results$ks_pval_raw_median, method = correction_method)
    corrected_q25 <- p.adjust(combined_results$ks_pval_raw_q25, method = correction_method)
    corrected_q75 <- p.adjust(combined_results$ks_pval_raw_q75, method = correction_method)
    
    # Add corrected p-values to results
    combined_results$ks_pval_corrected_median <- corrected_median
    combined_results$ks_pval_corrected_q25 <- corrected_q25
    combined_results$ks_pval_corrected_q75 <- corrected_q75
    combined_results$ks_correction_method <- correction_method
    
    return(combined_results)
  } else {
    return(data.frame())
  }
}

# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

# Interactive testing setup
if (interactive()) {
  source("src/utils.r")
  argvs$np <- "LOP_L"
  argvs$syn_type <- "pre"
  argvs$use_preset <- "temporal_known"
  argvs$ann <- "data/visual_neurons_anno.csv"
  argvs$meta <- "data/viz_meta.csv"
  argvs$preset <- "data/viz_preset.csv"
  argvs$sparse_limit <- 100L
  argvs$n_quantiles <- 100L
  argvs$n_bootstrap <- 1000L
  argvs$conf_int <- 95.0
  argvs$correction_method <- "bonferroni"
  argvs$ks_subsample_size <- 1000L
  argvs$ks_n_iterations <- 1000L
  argvs$ks_correction_method <- "fdr"
  argvs$cppsrc <- "./src/stats/bin/depth_stats.cpp"
  syn_path <- file.path("int/idv_mat/", paste0(argvs$np, "_rotated.csv.gz"))
} else {
  argvs$sparse_limit <- as.integer(argvs$sparse_limit)
  argvs$n_quantiles <- if (is.null(argvs$n_quantiles)) 100L else as.integer(argvs$n_quantiles)
  argvs$n_bootstrap <- if (is.null(argvs$n_bootstrap)) 1000L else as.integer(argvs$n_bootstrap)
  argvs$conf_int <- if (is.null(argvs$conf_int)) 95.0 else as.numeric(argvs$conf_int)
  argvs$correction_method <- if (is.null(argvs$correction_method)) "bonferroni" else argvs$correction_method
  argvs$ks_subsample_size <- if (is.null(argvs$ks_subsample_size)) 1000L else as.integer(argvs$ks_subsample_size)
  argvs$ks_n_iterations <- if (is.null(argvs$ks_n_iterations)) 1000L else as.integer(argvs$ks_n_iterations)
  argvs$ks_correction_method <- if (is.null(argvs$ks_correction_method)) "fdr" else argvs$ks_correction_method
  syn_path <- argvs$synf

  # Robust file path resolution for C++ source
  if (is.null(argvs$cppsrc)) {
    # Try multiple potential paths
    candidates <- c(
      "./depth_stats.cpp",
      "src/stats/bin/depth_stats.cpp"
    )
    found <- candidates[file.exists(candidates)][1]
    if (is.na(found)) {
      stop("Cannot find depth_stats.cpp in expected locations: ",
           paste(candidates, collapse=", "))
    }
    argvs$cppsrc <- found
  }
}

# Verify C++ source file exists before sourcing
if (!file.exists(argvs$cppsrc)) {
  stop(sprintf("C++ source file not found: %s", argvs$cppsrc))
}

Rcpp::sourceCpp(argvs$cppsrc)

# Load plot metadata and presets with validation
plot_meta_all <- validate_config_file(
  argvs$meta,
  c(
    "neuropil", "x_axis", "y_axis", "axis_1_func", "axis_2_func", "min1",
    "max1", "min2", "max2", "zid", "zinvert", "minz", "maxz", "outlayout",
    "outr1", "outr2", "outd1", "outd2"
  )
)
plot_meta <- plot_meta_all[neuropil == argvs$np]
if (nrow(plot_meta) == 0) stop(sprintf("No metadata found for neuropil: %s", argvs$np))

preset_all <- validate_config_file(
  argvs$preset,
  c(
    "preset", "palette", "color_guide", "filter_func", "color_by",
    "notch_split", "do_highlight", "hl_type", "hl_col", "hl_val",
    "known_only"
  )
)
preset <- preset_all[preset == argvs$use_preset]
if (nrow(preset) == 0) stop(sprintf("No preset found: %s", argvs$use_preset))

# Validate that required functions exist
validate_functions(preset, plot_meta)

# Load annotations and coordinates
opc_anno <- validate_config_file(argvs$ann, c("cell_type"))
np_coord <- fread(syn_path)

# Apply same filtering as visualization
filter_func <- get(preset$filter_func)
np_coord <- filter_type(
  np_coord, syn_type = argvs$syn_type, sparse_limit = argvs$sparse_limit
)
np_coord <- filter_func(np_coord, opc_anno, syn_type = argvs$syn_type)

# Keep known-only statistics consistent with the corresponding visualizations.
if (preset$known_only) {
  np_coord <- np_coord[newly_ann == "N"]
}

# Handle highlighting if specified
if (preset$do_highlight && preset$hl_type == "label") {
  np_coord[, highlight := get(preset$hl_col) == preset$hl_val]
}
if (preset$do_highlight && preset$hl_type == "path") {
  if (file.exists(preset$hl_val)) {
    hl_labels <- readLines(preset$hl_val)
    np_coord[, highlight := get(preset$hl_col) %in% hl_labels]
  }
}

# Split by Notch if required
if (preset$notch_split) {
  np_coord$notch_ntype <- paste(np_coord$Notch, np_coord$ntype, sep = "_")
  np_coord_list <- split(np_coord, np_coord$notch_ntype, drop = TRUE)
} else {
  np_coord_list <- list(all = np_coord)
}

# Calculate Quantile Wasserstein distance for each split
all_test_results <- list()

# If grouping by cell_type, we need to decide whether it is pre- or post.
if (preset$color_by == "cell_type") {
  preset$color_by <- paste(argvs$syn_type, "type", sep = "_")
}

for (i in names(np_coord_list)) {
  current_data <- np_coord_list[[i]]
  
  # Skip if insufficient data
  if (nrow(current_data) <= argvs$sparse_limit) {
    next
  }
  
  # Prepare data for testing
  if (preset$do_highlight) {
    test_data <- current_data[highlight == TRUE, .(
      group = get(preset$color_by),
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  } else {
    test_data <- current_data[, .(
      group = get(preset$color_by),
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  }
  
  # Apply z-axis inversion if needed
  if (plot_meta$zinvert) {
    test_data$depth <- test_data$depth * -1
  }
  
  # Compute Quantile Wasserstein distances
  wasserstein_results <- quantile_wasserstein_dist(
    test_data, 
    sparse_limit = argvs$sparse_limit,
    n_quantiles = argvs$n_quantiles,
    n_bootstrap = argvs$n_bootstrap,
    conf_int = argvs$conf_int
  )
  
  if (nrow(wasserstein_results) == 0) {
    next
  }
  
  # Perform subsampled KS tests
  ks_results <- perform_subsampled_ks_tests(
    test_data,
    sparse_limit = argvs$sparse_limit,
    subsample_size = argvs$ks_subsample_size,
    n_iterations = argvs$ks_n_iterations,
    correction_method = argvs$ks_correction_method
  )
  
  # Combine results
  if (nrow(ks_results) > 0) {
    # Merge Wasserstein and KS results by group pairs
    combined_results <- merge(wasserstein_results, ks_results, 
                              by = c("group1", "group2"), all.x = TRUE)
  } else {
    # If no KS results, add empty KS columns
    combined_results <- wasserstein_results
    combined_results$ks_pval_raw_median <- NA
    combined_results$ks_pval_raw_q25 <- NA
    combined_results$ks_pval_raw_q75 <- NA
    combined_results$ks_pval_corrected_median <- NA
    combined_results$ks_pval_corrected_q25 <- NA
    combined_results$ks_pval_corrected_q75 <- NA
    combined_results$ks_correction_method <- NA
    combined_results$ks_subsample_size <- NA
    combined_results$ks_n_iterations <- NA
  }
  
  # Add metadata
  combined_results$split <- i
  combined_results$neuropil <- argvs$np
  combined_results$syn_type <- argvs$syn_type
  combined_results$preset <- argvs$use_preset
  all_test_results[[i]] <- combined_results
}

# Combine all results
# Generate output filename
out_prefix <- paste(
  argvs$np, argvs$syn_type, argvs$use_preset,
  sep = "_"
)
output_file <- paste0(out_prefix, "_depth_stats.csv")

if (length(all_test_results) > 0) {
  final_results <- do.call(rbind, all_test_results)
  
  # Write results
  cat("Depth statistics (Wasserstein + KS) results written to:", output_file, "\n")
  cat("Number of comparisons:", nrow(final_results), "\n")
  cat("Confidence interval:", unique(final_results$conf_level), "%\n")
  
  # Print effect size summary
  median_effects <- final_results$wasserstein_median
  cat("Wasserstein effect size summary:\n")
  cat("  Min:", round(min(median_effects, na.rm = TRUE), 4), "\n")
  cat("  Max:", round(max(median_effects, na.rm = TRUE), 4), "\n")
  cat("  Mean:", round(mean(median_effects, na.rm = TRUE), 4), "\n")
  
  # Print KS p-value summary if available
  if ("ks_pval_raw_median" %in% colnames(final_results)) {
    ks_pvals <- final_results$ks_pval_raw_median
    valid_ks <- ks_pvals[!is.na(ks_pvals)]
    if (length(valid_ks) > 0) {
      cat("KS raw p-value summary:\n")
      cat("  Min:", round(min(valid_ks), 4), "\n")
      cat("  Max:", round(max(valid_ks), 4), "\n")
      cat("  Mean:", round(mean(valid_ks), 4), "\n")
      cat("KS subsample size:", unique(final_results$ks_subsample_size)[1], "\n")
      cat("KS iterations:", unique(final_results$ks_n_iterations)[1], "\n")
    }
  }
  
} else {
  final_results <- data.frame(
    group1 = logical(0),
    group2 = logical(0),
    wasserstein_median = logical(0),
    wasserstein_lower = logical(0),
    wasserstein_upper = logical(0),
    conf_level = logical(0),
    ks_pval_raw_median = logical(0),
    ks_pval_raw_q25 = logical(0),
    ks_pval_raw_q75 = logical(0),
    ks_pval_corrected_median = logical(0),
    ks_pval_corrected_q25 = logical(0),
    ks_pval_corrected_q75 = logical(0),
    ks_correction_method = logical(0),
    ks_subsample_size = logical(0),
    ks_n_iterations = logical(0),
    split = logical(0),
    neuropil = logical(0),
    syn_type = logical(0),
    preset = logical(0)
  )
  cat("No valid comparisons could be performed\n")
}

fwrite(final_results, output_file)
