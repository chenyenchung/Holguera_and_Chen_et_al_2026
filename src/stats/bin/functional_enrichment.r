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
suppressPackageStartupMessages(library(ggpattern))
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
  test_type, subsystem, n_total, n_subsys, test_used, pval
) {
  pval_collector[[length(pval_collector) + 1]] <<- list(
    test_type = test_type,
    subsystem = subsystem,
    n_total = n_total,
    n_subsys = n_subsys,
    test_used = test_used,
    pval = pval
  )
}

# Write all p-values with corrections to Excel
write_corrected_pvals_excel <- function(filename) {
  if (length(pval_collector) == 0) return()
  
  # Convert to data frame
  pval_df <- do.call(rbind.data.frame, pval_collector)
  
  # Apply FDR correction
  pval_df$p_fdr <- p.adjust(pval_df$pval, method = "fdr")
  
  # Rename columns for clarity
  colnames(pval_df) <- c("Test_Type", "Functional_Subsystem", "N_Total", 
                         "N_Subsystem", "Test_Used", "P_value_raw", "P_value_FDR")
  
  # Split data by test type
  test_types <- unique(pval_df$Test_Type)
  
  # Create workbook
  wb <- createWorkbook()
  
  # Add each test type as a separate sheet
  for (test_type in test_types) {
    sheet_data <- pval_df[pval_df$Test_Type == test_type, ]
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

# ============================================================================
# MODULAR ANALYSIS FUNCTIONS
# ============================================================================

# Run all statistical tests for a subsystem
run_all_tests <- function(data, subsys, el_cut) {
  results <- list()
  
  # Temporal association test
  tbl_temporal <- table(data$temporal_label, data$is_subsystem)
  if (nrow(tbl_temporal) >= 2 && ncol(tbl_temporal) >= 2) {
    test_obj <- tryCatch(
      {fisher.test(tbl_temporal)}, 
      error = function(e) {list(p.value = NA)}
    )
    results$temporal <- list(
      test = test_obj,
      n_total = sum(tbl_temporal),
      n_subsys = sum(tbl_temporal[, "TRUE"])
    )
    record_pval(
      "Temporal", subsys, results$temporal$n_total, 
      results$temporal$n_subsys, "fisher", test_obj$p.value
    )
  }
  
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

# Prepare data for temporal plot with hatched patterns
prepare_temporal_plot_data <- function(data, subsys) {
  stacked_data <- data %>%
    filter(Notch %in% c("Notch Off", "Notch On")) %>%
    group_by(temporal_label) %>%
    summarise(
      not_in_subsystem_notch_off = sum(!is_subsystem & Notch == "Notch Off") / n(),
      not_in_subsystem_notch_on = sum(!is_subsystem & Notch == "Notch On") / n(),
      in_subsystem_notch_off = sum(is_subsystem & Notch == "Notch Off") / n(),
      in_subsystem_notch_on = sum(is_subsystem & Notch == "Notch On") / n(),
      .groups = 'drop'
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        not_in_subsystem_notch_off, not_in_subsystem_notch_on,
        in_subsystem_notch_off, in_subsystem_notch_on
      ),
      names_to = "category",
      values_to = "prop"
    ) %>%
    mutate(
      category = factor(
        category,
        levels = c(
          "not_in_subsystem_notch_off",
          "not_in_subsystem_notch_on",
          "in_subsystem_notch_off",
          "in_subsystem_notch_on"
        )
      ),
      fill_color = factor(
        case_when(
          grepl("not_in_subsystem", category) ~ "Not in subsystem",
          TRUE ~ subsys
        ),
        levels = c("Not in subsystem", subsys)
      ),
      pattern_type = case_when(
        grepl("notch_off", category) ~ "none",
        grepl("notch_on", category) ~ "stripe"
      ),
      temporal_label = factor(
        temporal_label,
        levels = unique(data$temporal_label)
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

# Create temporal plot with hatched patterns
create_temporal_hatched_plot <- function(
    plot_data, subsys, test_results, pval_df = NULL) {
  # Get temporal test p-value
  p_value <- ifelse(
    is.null(test_results$temporal$test$p.value), 
    NA, 
    test_results$temporal$test$p.value
  )
  
  # Get FDR-corrected p-value if available
  p_fdr <- p_value
  if (!is.null(pval_df)) {
    # Find the FDR-corrected p-value for this specific test
    idx <- which(pval_df$test_type == "Temporal" & pval_df$subsystem == subsys)
    if (length(idx) > 0) {
      p_fdr <- pval_df$p_fdr[idx[1]]
    }
  }
  
  p <- ggplot(plot_data, aes(x = temporal_label, y = prop)) +
    geom_col_pattern(
      aes(fill = fill_color, pattern = pattern_type),
      position = position_stack(),
      pattern_fill = "grey90",
      pattern_colour = "black",
      pattern_density = 0.3,
      pattern_spacing = 0.1,
      pattern_size = 0.1,
      pattern_key_scale_factor = 0.2
    ) +
    labs(
      title = subsys,
      subtitle = paste0("Fisher's Exact Test (FDR) p = ", 
                       round(p_fdr, 3)),
      x = "Temporal Origin",
      y = "Proportion of neurons",
      fill = paste("Involved in\n", subsys),
      pattern = "Notch Status"
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_fill_manual(
      values = setNames(
        c("grey90", scale_fill_subsystem()$palette(0)[[subsys]]),
        c("Not in subsystem", subsys)
      ),
      labels = c("No", "Yes")
    ) +
    scale_pattern_manual(
      values = c("none" = "none", "stripe" = "stripe"),
      labels = c("Off", "On")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      legend.position = "right",
      legend.box = "vertical",
      axis.title.x = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.subtitle = element_text(size = 10, color = "grey50"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    guides(
      pattern = guide_legend(
        override.aes = list(fill = scale_fill_subsystem()$palette(0)[[subsys]]),
        order = 1
        ),
      fill = guide_legend(override.aes = list(pattern = "none"),
                          order = 2)
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
hatched_plots <- list()
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
  
  # Prepare data for temporal plot with hatched patterns
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
pval_df <- NULL
if (length(pval_collector) > 0) {
  pval_df <- do.call(rbind.data.frame, pval_collector)
  pval_df$p_fdr <- p.adjust(pval_df$pval, method = "fdr")
}

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
  
  # Create temporal plot with FDR-corrected p-values
  if (!is.null(plot_data$temporal)) {
    hatched_plots[[subsys]] <- create_temporal_hatched_plot(
      plot_data$temporal, subsys, test_results, pval_df
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

# Save temporal plots with hatched patterns
if (length(hatched_plots) > 0) {
  outfn <- "functional_subsystems_temporal_by_notch_stacked.pdf"
  wrap_plots(hatched_plots, ncol = 2)
  ggsave(
    outfn, 
    width = 16, 
    height = 4 * ceiling(length(hatched_plots)/2)
  )
  message(sprintf("Saved temporal hatched plots to: %s", outfn))
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

pval_summary <- write_corrected_pvals_excel(pval_file)

# Print summary of significant results
message("\n=== SUMMARY OF RESULTS (FDR-corrected) ===")
if (!is.null(pval_summary)) {
  sig_results <- pval_summary[pval_summary$p_fdr < 0.05, ]
  if (nrow(sig_results) > 0) {
    message("Significant associations after FDR correction:")
    for (i in 1:nrow(sig_results)) {
      message(sprintf("  %s - %s: p(FDR) = %.4f", 
                     sig_results$test_type[i],
                     sig_results$subsystem[i], 
                     sig_results$p_fdr[i]))
    }
  } else {
    message("No significant associations after FDR correction")
  }
}

message("\nFunctional enrichment analysis complete!")
message("Results saved to ./")
message("Corrected p-values saved to: ", pval_file)
