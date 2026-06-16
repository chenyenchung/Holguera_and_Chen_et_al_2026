#!/usr/bin/env Rscript
renv::load("/scratch/ycc520/flyem")
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

# ============================================================================
# UTILITY FUNCTIONS FOR DEPTH STATS VISUALIZATION
# ============================================================================

#' Convert p-values to significance stars
#' @param p Numeric p-value
#' @return Character string with significance annotation
p2star <- function(p) {
  if (is.na(p)) return("")
  if (p > 0.05) {
    return("ns")
  } else if (p <= 0.0001) {
    return("****")
  } else if (p <= 0.001) {
    return("***")
  } else if (p <= 0.01) {
    return("**")
  } else {
    return("*")
  }
}

#' Convert pairwise comparison data to matrix format
#' @param data Data frame with group1, group2, and value columns
#' @param value_col Name of the column containing values
#' @return Symmetric matrix
pairwise_to_matrix <- function(data, value_col, diag_zero = TRUE) {
  # Get all unique groups
  all_groups <- unique(c(data$group1, data$group2))
  
  # Create empty matrix
  mat <- matrix(NA, nrow = length(all_groups), ncol = length(all_groups))
  rownames(mat) <- all_groups
  colnames(mat) <- all_groups
  
  if (diag_zero) {
    # Fill diagonal with zeros (self-comparison)
    diag(mat) <- 0
  } else {
    diag(mat) <- NA
  }
  
  # Fill matrix with values
  for (i in 1:nrow(data)) {
    g1 <- data$group1[i]
    g2 <- data$group2[i]
    val <- data[[value_col]][i]
    
    mat[g1, g2] <- val
    mat[g2, g1] <- val  # Make symmetric
  }
  
  return(mat)
}

#' Create matrix heatmap
#' @param mat Matrix to plot
#' @param title Plot title
#' @param color_scale Color scale function
#' @param show_values Whether to show values in cells
#' @param value_labels Optional matrix of labels to show in cells
#' @return ggplot object
create_matrix_heatmap <- function(mat, title, color_scale = NULL, 
                                  show_values = FALSE, 
                                  value_labels = NULL) {
  
  # Convert matrix to long format
  mat_long <- melt(mat, varnames = c("group1", "group2"), value.name = "value")
  
  # Create base plot
  p <- ggplot(mat_long, aes(x = group1, y = group2, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    labs(title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 0, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )
  
  # Add color scale
  if (!is.null(color_scale)) {
    p <- p + color_scale
  }
  
  # Add value labels if requested
  if (show_values && !is.null(value_labels)) {
    # Convert value_labels matrix to long format
    labels_long <- melt(value_labels, varnames = c("group1", "group2"), value.name = "label")
    
    p <- p +
      geom_text(data = labels_long, aes(x = group1, y = group2, label = label), 
                color = "black", size = 5, inherit.aes = FALSE)
  }
  
  return(p)
}

# ============================================================================
# PLOTTING FUNCTIONS FOR DEPTH STATS
# ============================================================================

#' Plot 95% CI of Wasserstein distance
#' @param data Data frame with Wasserstein test results
#' @return ggplot object
plot_wasserstein_ci <- function(data) {
  # Create comparison labels
  data$comparison <- paste(data$group1, "vs", data$group2, sep = "\n")
  
  p <- ggplot(data, aes(x = comparison, y = wasserstein_median)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = wasserstein_lower, ymax = wasserstein_upper),
                  width = 0.2) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(
      title = "95% CIs of Synaptic Depth Distribution Difference",
      y = "Wasserstein Distance"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  return(p)
}

#' Plot paired heatmap of Wasserstein distance
#' @param data Data frame with Wasserstein test results
#' @return ggplot object
plot_wasserstein_heatmap <- function(data) {
  # Convert to matrix
  w_matrix <- pairwise_to_matrix(data, "wasserstein_median")
  
  # Create heatmap
  p <- create_matrix_heatmap(
    w_matrix,
    title = "Pairwise Synaptic Depth Distribution Difference",
    color_scale = scale_fill_gradient(
      low = "white", high = "blue",
      name = "Wasserstein\nDistance"
    )
  )
  
  return(p)
}

#' Plot point-and-whisker plot of K-S p-values
#' @param data Data frame with K-S test results
#' @return ggplot object
plot_ks_boxplot <- function(data) {
  # Create comparison labels
  data$comparison <- paste(data$group1, "vs", data$group2, sep = "\n")
  
  # Prepare data for boxplot-style visualization
  plot_data <- data %>%
    select(comparison, ks_pval_raw_median, ks_pval_raw_q25, ks_pval_raw_q75) %>%
    filter(!is.na(ks_pval_raw_median))
  
  p <- ggplot(plot_data, aes(x = comparison, y = ks_pval_raw_median)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ks_pval_raw_q25, ymax = ks_pval_raw_q75),
                  width = 0.2) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = "Subsampled Kolmogorov-Smirnov P-values",
      y = "P-value (median with quartiles)"
    ) +
    scale_y_continuous(transform = c("log10", "reverse")) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  return(p)
}

#' Plot significance heatmap of K-S p-values
#' @param data Data frame with K-S test results
#' @return ggplot object
plot_ks_significance_heatmap <- function(data) {
  # Filter out missing p-values
  data_clean <- data %>% filter(!is.na(ks_pval_corrected_median))
  
  if (nrow(data_clean) == 0) {
    warning("No valid corrected p-values found")
    return(NULL)
  }
  
  # Convert to matrix
  p_matrix <- pairwise_to_matrix(data_clean, "ks_pval_corrected_median", diag_zero = FALSE)
  
  # Handle log transformation with edge case p = 0 and NA values
  log_p_matrix <- p_matrix
  
  # Preserve NA values
  na_mask <- is.na(p_matrix)
  
  # Handle edge case where p = 0
  zero_mask <- !na_mask & p_matrix == 0
  positive_mask <- !na_mask & p_matrix > 0
  
  # Calculate minimum positive value, with fallback
  min_positive <- if (any(positive_mask)) {
    min(p_matrix[positive_mask], na.rm = TRUE) / 10
  } else {
    1e-10  # Fallback value when no positive p-values exist
  }
  
  # Apply transformations
  log_p_matrix[zero_mask] <- -log10(min_positive)
  log_p_matrix[positive_mask] <- -log10(p_matrix[positive_mask])
  log_p_matrix[na_mask] <- NA
  
  # Create significance labels matrix
  sig_matrix <- matrix(NA, nrow = nrow(p_matrix), ncol = ncol(p_matrix))
  rownames(sig_matrix) <- rownames(p_matrix)
  colnames(sig_matrix) <- colnames(p_matrix)
  
  for (i in 1:nrow(p_matrix)) {
    for (j in 1:ncol(p_matrix)) {
      if (!is.na(p_matrix[i, j])) {
        sig_matrix[i, j] <- p2star(p_matrix[i, j])
      }
    }
  }
  
  # Create heatmap
  p <- create_matrix_heatmap(
    log_p_matrix,
    title = "K-S P-values (Corrected)",
    color_scale = scale_fill_gradient(
      low = "white", high = "red",
      name = "Adjusted\np-value",
      labels = function(x) paste0("1e-", x)
    ),
    show_values = TRUE,
    value_labels = sig_matrix
  )
  
  return(p)
}

# ============================================================================
# MAIN VISUALIZATION FUNCTION
# ============================================================================

#' Generate comprehensive depth statistics visualization report
#' @param test_results Data frame with depth stats results 
#' @param output_prefix Output file prefix
#' @return List of generated file paths
generate_depth_stats_report <- function(test_results, output_prefix) {
  
  generated_files <- c()
  
  # Calculate dynamic dimensions
  n_groups <- length(unique(c(test_results$group1, test_results$group2)))
  n_comparisons <- nrow(test_results)
  
  # Dynamic sizing for different plot types
  point_whisker_width <- max(6, n_comparisons * 1 + 1)
  point_whisker_height <- 6
  heatmap_size <- max(6, n_groups * 0.6 + 1)
  
  # Generate 95% CI plot for Wasserstein distances
  cat("Generating Wasserstein CI plot...\n")
  ci_plot <- plot_wasserstein_ci(test_results)
  ci_file <- paste0(output_prefix, "_wasserstein_ci.pdf")
  ggsave(plot = ci_plot, filename = ci_file, width = point_whisker_width, height = point_whisker_height, dpi = 300, limitsize = FALSE)
  generated_files <- c(generated_files, ci_file)
  
  # Generate Wasserstein distance heatmap
  cat("Generating Wasserstein heatmap...\n")
  w_heatmap <- plot_wasserstein_heatmap(test_results)
  w_file <- paste0(output_prefix, "_wasserstein_heatmap.pdf")
  ggsave(plot = w_heatmap, filename = w_file, width = heatmap_size + 1, height = heatmap_size, dpi = 300, limitsize = FALSE)
  generated_files <- c(generated_files, w_file)
  
  # Generate K-S plots if data is available
  if ("ks_pval_raw_median" %in% colnames(test_results)) {
    cat("Generating K-S boxplot...\n")
    ks_boxplot <- plot_ks_boxplot(test_results)
    ks_file <- paste0(output_prefix, "_ks_boxplot.pdf")
    ggsave(plot = ks_boxplot, filename = ks_file, width = point_whisker_width, height = point_whisker_height, dpi = 300, limitsize = FALSE)
    generated_files <- c(generated_files, ks_file)
    
    if ("ks_pval_corrected_median" %in% colnames(test_results)) {
      cat("Generating K-S significance heatmap...\n")
      ks_heatmap <- plot_ks_significance_heatmap(test_results)
      if (!is.null(ks_heatmap)) {
        ks_hm_file <- paste0(output_prefix, "_ks_significance_heatmap.pdf")
        ggsave(plot = ks_heatmap, filename = ks_hm_file, width = heatmap_size + 1, height = heatmap_size, dpi = 300, limitsize = FALSE)
        generated_files <- c(generated_files, ks_hm_file)
      }
    }
  }
  
  return(generated_files)
}

# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  argvs$input_file <- list.files(pattern = ".*\\.csv")[[1]]
  argvs$output_prefix <- "test_depth_viz"
  source("./src/utils.r")
}

# Check if running as standalone script
if (!is.null(argvs$input_file)) {
  
  # Load test results
  if (!file.exists(argvs$input_file)) {
    stop("Input file not found: ", argvs$input_file)
  }
  
  test_results <- fread(argvs$input_file)
  
  # Validate required columns
  required_cols <- c("group1", "group2", "wasserstein_median", 
                     "wasserstein_lower", "wasserstein_upper")
  missing_cols <- setdiff(required_cols, colnames(test_results))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Set defaults
  output_prefix <- if (is.null(argvs$output_prefix)) {
    gsub("_depth_stats\\.csv$", "", basename(argvs$input_file))
  } else {
    argvs$output_prefix
  }
  
  
  # Check if split column exists and split data accordingly
  if ("split" %in% colnames(test_results)) {
    cat("Split column detected. Splitting data by:", paste(unique(test_results$split), collapse = ", "), "\n")
    test_results_list <- split(test_results, test_results$split, drop = TRUE)
  } else {
    test_results_list <- list(all = test_results)
  }
  
  # Generate plots for each split
  cat("Loading data from:", argvs$input_file, "\n")
  cat("Data loaded successfully with", nrow(test_results), "comparisons\n")
  
  all_generated_files <- c()
  
  for (split_name in names(test_results_list)) {
    split_data <- test_results_list[[split_name]]
    
    # Create output prefix with split name
    if (split_name != "all") {
      split_output_prefix <- paste(output_prefix, split_name, sep = "_")
      cat("\nGenerating plots for split:", split_name, "\n")
    } else {
      split_output_prefix <- output_prefix
    }
    
    # Generate plots for this split
    generated_files <- generate_depth_stats_report(split_data, split_output_prefix)
    all_generated_files <- c(all_generated_files, generated_files)
  }
  
  cat("\nGenerated files:\n")
  for (file in all_generated_files) {
    cat("  ", file, "\n")
  }
  
  # Print summary statistics
  cat("\nSummary Statistics:\n")
  cat("Total comparisons:", nrow(test_results), "\n")
  if ("split" %in% colnames(test_results)) {
    cat("Number of splits:", length(test_results_list), "\n")
    for (split_name in names(test_results_list)) {
      cat("  ", split_name, ":", nrow(test_results_list[[split_name]]), "comparisons\n")
    }
  }
  
  if ("conf_level" %in% colnames(test_results)) {
    cat("Confidence level:", unique(test_results$conf_level), "%\n")
  }
  
  # Wasserstein distance summary
  effect_sizes <- test_results$wasserstein_median
  cat("Wasserstein distance range:", round(min(effect_sizes, na.rm = TRUE), 4), "to", 
      round(max(effect_sizes, na.rm = TRUE), 4), "\n")
  
  # K-S p-value summary if available
  if ("ks_pval_raw_median" %in% colnames(test_results)) {
    ks_pvals <- test_results$ks_pval_raw_median
    valid_ks <- ks_pvals[!is.na(ks_pvals)]
    if (length(valid_ks) > 0) {
      cat("K-S raw p-value range:", round(min(valid_ks), 4), "to", 
          round(max(valid_ks), 4), "\n")
    }
  }
  
  cat("\nVisualization complete!\n")
}
