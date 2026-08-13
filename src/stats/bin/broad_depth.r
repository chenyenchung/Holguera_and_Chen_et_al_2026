#!/usr/bin/env Rscript
renv::load("/scratch/ycc520/flyem")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Rcpp))

perform_broad_depth_analysis <- function(data, ref_sup_mapping, ref_deep_mapping,
                                       types_of_interest, coefficient = 0.5,
                                       n_bootstrap = 1000, conf_int = 95,
                                       seed = 1L, syn_type = "pre",
                                       ref_superficial, ref_deep,
                                       data_superficial, data_deep) {

  interest_data <- data[group %in% types_of_interest]

  # Prepare interest group neuron mapping
  # For interest group: only one role (pre or post depending on syn_type)
  interest_col <- if (syn_type == "pre") "pre_root_id" else "post_root_id"
  depth_col <- "depth"  # Already set to correct depth in test_data

  interest_data_prep <- interest_data[, .(
    neuron_id = get(interest_col),
    depth = get(depth_col)
  )]
  setorder(interest_data_prep, neuron_id)
  interest_neuron_info <- interest_data_prep[, .N, by = neuron_id]
  interest_neuron_info[, start_idx := c(0, cumsum(N[-.N]))]

  interest_mapping <- list(
    depths = interest_data_prep$depth,
    neuron_starts = interest_neuron_info$start_idx,
    neuron_counts = interest_neuron_info$N
  )

  # Call C++ neuron-level bootstrap function
  bootstrap_result <- perform_broad_depth_bootstrap_neuron_level(
    ref_sup_depths = ref_sup_mapping$depths,
    ref_sup_neuron_starts = as.integer(ref_sup_mapping$neuron_starts),
    ref_sup_neuron_counts = as.integer(ref_sup_mapping$neuron_counts),

    ref_deep_depths = ref_deep_mapping$depths,
    ref_deep_neuron_starts = as.integer(ref_deep_mapping$neuron_starts),
    ref_deep_neuron_counts = as.integer(ref_deep_mapping$neuron_counts),

    interest_depths = interest_mapping$depths,
    interest_neuron_starts = as.integer(interest_mapping$neuron_starts),
    interest_neuron_counts = as.integer(interest_mapping$neuron_counts),

    coefficient = coefficient,
    n_bootstrap = n_bootstrap,
    conf_int = conf_int,
    seed = if (is.null(seed)) NULL else as.integer(seed)
  )
  
  # Determine observed direction based on bootstrap results
  if (abs(bootstrap_result$observed_distance_diff) > bootstrap_result$observed_delta_thres) {
    observed_direction <- if (bootstrap_result$observed_distance_diff > 0) "superficial" else "deep"
  } else {
    observed_direction <- "neither"
  }
  
  result <- list(
    types_of_interest = paste(types_of_interest, collapse = ";"),
    ref_superficial = paste(ref_superficial, collapse = ";"),
    ref_deep = paste(ref_deep, collapse = ";"),
    n_types_interest = length(types_of_interest),
    n_ref_superficial = length(ref_superficial),
    n_ref_deep = length(ref_deep),
    n_samples_interest = nrow(interest_data),
    n_samples_superficial = nrow(data_superficial),
    n_samples_deep = nrow(data_deep),
    n_neurons_interest = bootstrap_result$n_neurons_interest,
    n_neurons_superficial = bootstrap_result$n_neurons_sup,
    n_neurons_deep = bootstrap_result$n_neurons_deep,
    observed_sup_d = bootstrap_result$observed_sup_d,
    observed_deep_d = bootstrap_result$observed_deep_d,
    observed_distance_diff = bootstrap_result$observed_distance_diff,
    observed_delta_thres_base = bootstrap_result$observed_delta_thres_base,
    observed_delta_thres = bootstrap_result$observed_delta_thres,
    observed_test_statistic = bootstrap_result$observed_test_statistic,
    observed_bias_ratio = bootstrap_result$observed_bias_ratio,
    bootstrap_distance_diff_median = bootstrap_result$bootstrap_distance_diff_median,
    bootstrap_distance_diff_lower = bootstrap_result$bootstrap_distance_diff_lower,
    bootstrap_distance_diff_upper = bootstrap_result$bootstrap_distance_diff_upper,
    bootstrap_delta_thres_base_median = bootstrap_result$bootstrap_delta_thres_base_median,
    bootstrap_delta_thres_base_lower = bootstrap_result$bootstrap_delta_thres_base_lower,
    bootstrap_delta_thres_base_upper = bootstrap_result$bootstrap_delta_thres_base_upper,
    bootstrap_test_stat_median = bootstrap_result$bootstrap_test_stat_median,
    bootstrap_test_stat_lower = bootstrap_result$bootstrap_test_stat_lower,
    bootstrap_test_stat_upper = bootstrap_result$bootstrap_test_stat_upper,
    bootstrap_bias_ratio_median = bootstrap_result$bootstrap_bias_ratio_median,
    bootstrap_bias_ratio_lower = bootstrap_result$bootstrap_bias_ratio_lower,
    bootstrap_bias_ratio_upper = bootstrap_result$bootstrap_bias_ratio_upper,
    coefficient = coefficient,
    p_value_exceeds_threshold = bootstrap_result$p_value_exceeds_threshold,
    p_value_bias_ratio = bootstrap_result$p_value_bias_ratio,
    direction = observed_direction,
    conf_level = conf_int,
    n_bootstrap = n_bootstrap
  )
  
  return(result)
}

# Derive a reproducible seed for each comparison so results do not depend on
# comparison order while still allowing the full run to be changed by one seed.
stable_seed <- function(base_seed, ...) {
  modulus <- 2147483647
  value <- as.numeric(base_seed) %% modulus
  key <- paste(..., sep = "|")

  for (byte in as.integer(charToRaw(key))) {
    value <- (value * 131 + byte) %% modulus
  }

  as.integer(max(1, value))
}

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}
if (file.exists("./src/utils.r")) {
  source("src/utils.r", chdir = FALSE)
}

if (interactive()) {
  source("src/utils.r")
  argvs$np <- "LOP_L"
  argvs$syn_type <- "pre"
  argvs$use_preset <- "type_putative"
  argvs$ann <- "data/visual_neurons_anno.csv"
  argvs$meta <- "data/viz_meta.csv"
  argvs$preset <- "data/viz_preset.csv"
  argvs$sparse_limit <- 100L
  argvs$coefficient <- 0.5
  argvs$n_bootstrap <- 1000L
  argvs$conf_int <- 95.0
  argvs$bootstrap_seed <- 1L
  argvs$cppsrc <- "src/stats/bin/broad_depth.cpp"
  syn_path <- file.path("int/idv_mat/", paste0(argvs$np, "_rotated.csv.gz"))
} else {
  argvs$sparse_limit <- as.integer(argvs$sparse_limit)
  argvs$coefficient <- if (is.null(argvs$coefficient)) 0.5 else as.numeric(argvs$coefficient)
  argvs$n_bootstrap <- if (is.null(argvs$n_bootstrap)) 1000L else as.integer(argvs$n_bootstrap)
  argvs$conf_int <- if (is.null(argvs$conf_int)) 95.0 else as.numeric(argvs$conf_int)
  argvs$bootstrap_seed <- if (is.null(argvs$bootstrap_seed)) 1L else as.integer(argvs$bootstrap_seed)
  syn_path <- argvs$synf
  if (is.null(argvs$preset)) argvs$preset <- "./viz_preset.csv"
  if (is.null(argvs$use_preset)) stop("use_preset parameter is required")

  # Robust file path resolution for C++ source
  if (is.null(argvs$cppsrc)) {
    # Try multiple potential paths
    candidates <- c(
      "./broad_depth.cpp",
      "src/stats/bin/broad_depth.cpp"
    )
    found <- candidates[file.exists(candidates)][1]
    if (is.na(found)) {
      stop("Cannot find broad_depth.cpp in expected locations: ",
           paste(candidates, collapse=", "))
    }
    argvs$cppsrc <- found
  }
}

if (length(argvs$bootstrap_seed) != 1 || is.na(argvs$bootstrap_seed)) {
  stop("bootstrap_seed must be a single integer")
}
cat("Broad-depth bootstrap base seed:", argvs$bootstrap_seed, "\n")

# Verify C++ source file exists before sourcing
if (!file.exists(argvs$cppsrc)) {
  stop(sprintf("C++ source file not found: %s", argvs$cppsrc))
}

# Compile C++ source
Rcpp::sourceCpp(argvs$cppsrc)

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

opc_anno <- validate_config_file(argvs$ann, c("cell_type"))
np_coord <- fread(syn_path)

# Store unfiltered copy for reference group calculation
np_coord_unfiltered <- copy(np_coord)

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

# Get reference depths
available_groups <- unique(
  c(np_coord_unfiltered$pre_type, np_coord_unfiltered$post_type)
)

# Load reference groups from configuration file
ref_groups_file <- if (!is.null(argvs$ref_groups)) {
  argvs$ref_groups
} else if (file.exists("data/reference_groups.csv")) {
  "data/reference_groups.csv"
} else if (file.exists("../data/reference_groups.csv")) {
  "../data/reference_groups.csv"
} else {
  NULL
}

if (!is.null(ref_groups_file) && file.exists(ref_groups_file)) {
  # Load reference groups from CSV
  ref_groups_config <- fread(ref_groups_file)
  ref_config <- ref_groups_config[neuropil == argvs$np]

  if (nrow(ref_config) > 0) {
    # Apply pattern based on type
    if (ref_config$pattern_type == "regex") {
      ref_superficial <- grep(ref_config$ref_superficial, available_groups, value = TRUE)
      ref_deep <- grep(ref_config$ref_deep, available_groups, value = TRUE)
    } else if (ref_config$pattern_type == "exact") {
      ref_superficial <- intersect(
        strsplit(ref_config$ref_superficial, ";")[[1]],
        available_groups
      )
      ref_deep <- intersect(
        strsplit(ref_config$ref_deep, ";")[[1]],
        available_groups
      )
    } else {
      warning("Unknown pattern_type in reference_groups.csv: ", ref_config$pattern_type)
      ref_superficial <- character(0)
      ref_deep <- character(0)
    }
  } else {
    warning("No reference groups found for neuropil: ", argvs$np)
    ref_superficial <- character(0)
    ref_deep <- character(0)
  }
} else {
  # Fallback to hardcoded reference groups with deprecation warning
  warning("Reference groups configuration file not found. Using hardcoded values. ",
          "Please create data/reference_groups.csv for future analyses.")

  if (grepl("^ME_", argvs$np)) {
    ref_superficial <- grep("^Dm", available_groups, value = TRUE)
    ref_deep <- grep("^Pm", available_groups, value = TRUE)
  } else if (grepl("^LOP_", argvs$np)) {
    ref_superficial <- grep("T(4|5)(a|b)$", available_groups, value = TRUE)
    ref_deep <- grep("T(4|5)(c|d)$", available_groups, value = TRUE)
  } else if (grepl("^LO_", argvs$np)) {
    ref_superficial <- c("Tm1", "Tm2", "Tm4")
    ref_deep <- c("Tm20", "Tm5a", "Tm5b", "Tm5c", "Tm5d", "Tm5e", "Tm5f")
  } else {
    # For other neuropils, use a generic approach or skip analysis
    ref_superficial <- character(0)
    ref_deep <- character(0)
  }
}

# Log reference group composition for validation
if (length(ref_superficial) > 0 || length(ref_deep) > 0) {
  cat(sprintf("\n=== Reference Groups for %s (syn_type=%s) ===\n", argvs$np, argvs$syn_type))
  cat(sprintf("Superficial (%d types): %s\n",
              length(ref_superficial),
              paste(sort(ref_superficial), collapse=", ")))
  cat(sprintf("Deep (%d types): %s\n",
              length(ref_deep),
              paste(sort(ref_deep), collapse=", ")))
  cat("==========================================\n\n")
}

# Helper function to create neuron-to-synapse mapping for neuron-level bootstrap
prepare_neuron_depth_mapping <- function(data, neuron_types, syn_type) {
  # For reference groups: combine pre & post neurons
  # Create "pre_<root_id>" and "post_<root_id>" to distinguish same neuron in different roles

  # Collect all neurons and their synapses' depths
  pre_part <- data[pre_type %in% neuron_types, .(
    neuron_id = paste0("pre_", pre_root_id),
    depth = pre_rz
  )]
  post_part <- data[post_type %in% neuron_types, .(
    neuron_id = paste0("post_", post_root_id),
    depth = post_rz
  )]
  all_data <- rbind(pre_part, post_part)

  # Sort by neuron_id to group synapses from same neuron
  setorder(all_data, neuron_id)

  # Create neuron grouping info: for each unique neuron, where its synapses start and how many
  neuron_info <- all_data[, .N, by = neuron_id]
  neuron_info[, start_idx := c(0, cumsum(N[-.N]))]  # 0-based indexing for C++

  return(list(
    depths = all_data$depth,
    neuron_starts = neuron_info$start_idx,
    neuron_counts = neuron_info$N
  ))
}

# Use unfiltered data for reference groups
# Reference types are selected landmarks and should include ALL their connections
data_superficial <- np_coord_unfiltered[
  pre_type %in% ref_superficial | post_type %in% ref_superficial
]
data_deep <- np_coord_unfiltered[
  pre_type %in% ref_deep | post_type %in% ref_deep
]

# Apply z-axis inversion BEFORE creating neuron mappings
# This ensures reference depths match the coordinate system used for test data
if (plot_meta$zinvert) {
  data_superficial$pre_rz <- data_superficial$pre_rz * -1
  data_superficial$post_rz <- data_superficial$post_rz * -1
  data_deep$pre_rz <- data_deep$pre_rz * -1
  data_deep$post_rz <- data_deep$post_rz * -1
}

# Create neuron-level mappings for reference groups
ref_sup_mapping <- prepare_neuron_depth_mapping(data_superficial, ref_superficial, argvs$syn_type)
ref_deep_mapping <- prepare_neuron_depth_mapping(data_deep, ref_deep, argvs$syn_type)

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

# Perform broad depth analysis for each split
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
  
  # Prepare test data using preset grouping
  if (preset$do_highlight) {
    test_data <- current_data[highlight == TRUE, .(
      group = get(preset$color_by),
      pre_type,
      post_type,
      pre_root_id,
      post_root_id,
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  } else {
    test_data <- current_data[, .(
      group = get(preset$color_by),
      pre_type,
      post_type,
      pre_root_id,
      post_root_id,
      depth = get(paste(argvs$syn_type, "rz", sep = "_"))
    )]
  }
  
  # Apply z-axis inversion if needed
  if (plot_meta$zinvert) {
    test_data$depth <- test_data$depth * -1
  }
  
  # Get groups to test
  groups_to_test <- unique(as.character(test_data$group))
  groups_to_test <- groups_to_test[!is.na(groups_to_test)]
  if (preset$color_by == "spatial_origin") {
    groups_to_test <- setdiff(groups_to_test, "unknown")
  }
  
  # Perform bootstrap analysis for each group
  if (length(ref_superficial) > 0 && length(ref_deep) > 0 && length(groups_to_test) > 0) {
    split_results <- lapply(groups_to_test, function(type) {
      analysis_seed <- stable_seed(
        argvs$bootstrap_seed,
        argvs$np,
        argvs$syn_type,
        argvs$use_preset,
        i,
        type
      )
      out <- perform_broad_depth_analysis(
        test_data,
        ref_sup_mapping = ref_sup_mapping,
        ref_deep_mapping = ref_deep_mapping,
        types_of_interest = type,
        coefficient = argvs$coefficient,
        n_bootstrap = argvs$n_bootstrap,
        conf_int = argvs$conf_int,
        seed = analysis_seed,
        syn_type = argvs$syn_type,
        ref_superficial = ref_superficial,
        ref_deep = ref_deep,
        data_superficial = data_superficial,
        data_deep = data_deep
      )
      return(out)
    })
    
    # Combine results for this split
    if (length(split_results) > 1) {
      split_result_df <- do.call(rbind.data.frame, split_results)
      split_result_df$split <- i
      all_test_results[[i]] <- split_result_df
    } else if (length(split_results) == 1) {
      split_result_df <- as.data.frame(split_results[[1]])
      split_result_df$split <- i
      all_test_results[[i]] <- split_result_df
    } else {
      warning(paste("Insufficient reference groups or groups to test for bootstrap analysis in split:", i))
    }
  }
}

# Combine all results from all splits
if (length(all_test_results) > 0) {
  result <- do.call(rbind, all_test_results)
  
  # Add common metadata if not already present
  if (!"neuropil" %in% colnames(result)) {
    result$neuropil <- argvs$np
  }
  if (!"syn_type" %in% colnames(result)) {
    result$syn_type <- argvs$syn_type
  }
  if (!"preset" %in% colnames(result)) {
    result$preset <- argvs$use_preset
  }
  
  # Apply FDR correction for multiple testing, per-split
  # This ensures each split (e.g., Off_Interneurons vs On_Interneurons) is corrected independently
  if ("p_value_exceeds_threshold" %in% colnames(result)) {
    # Initialize FDR columns
    result$p_value_fdr <- NA_real_
    result$significant_fdr <- NA

    # Split results by the 'split' column to apply FDR correction within each split
    if ("split" %in% colnames(result)) {
      split_groups <- split(result, result$split)

      # Apply FDR correction to each split separately
      result_list <- lapply(split_groups, function(split_df) {
        p_values <- split_df$p_value_exceeds_threshold
        valid_indices <- !is.na(p_values)

        if (sum(valid_indices) > 1) {
          # Multiple valid p-values: apply FDR correction
          valid_p <- p_values[valid_indices]
          fdr_corrected <- p.adjust(valid_p, method = "fdr")

          # Fill in corrected values
          split_df$p_value_fdr[valid_indices] <- fdr_corrected
          split_df$significant_fdr[valid_indices] <- fdr_corrected < 0.05
        } else if (sum(valid_indices) == 1) {
          # Single valid p-value: FDR correction equals raw p-value
          split_df$p_value_fdr[valid_indices] <- p_values[valid_indices]
          split_df$significant_fdr[valid_indices] <- p_values[valid_indices] < 0.05
        }
        # If no valid p-values (sum == 0), leave as NA (already initialized)

        return(split_df)
      })

      # Recombine all splits
      result <- do.call(rbind, result_list)
    } else {
      # No split column - apply FDR to entire dataset
      p_values <- result$p_value_exceeds_threshold
      valid_indices <- !is.na(p_values)

      if (sum(valid_indices) > 1) {
        valid_p <- p_values[valid_indices]
        fdr_corrected <- p.adjust(valid_p, method = "fdr")
        result$p_value_fdr[valid_indices] <- fdr_corrected
        result$significant_fdr[valid_indices] <- fdr_corrected < 0.05
      } else if (sum(valid_indices) == 1) {
        result$p_value_fdr[valid_indices] <- p_values[valid_indices]
        result$significant_fdr[valid_indices] <- p_values[valid_indices] < 0.05
      }
    }
  }
} else {
  result <- NULL
}

out_prefix <- paste(
  argvs$np, argvs$syn_type, argvs$use_preset,
  sep = "_"
)
output_file <- paste0(out_prefix, "_broad_depth.csv")

if (!is.null(result) && nrow(result) > 0) {
  
  cat("Broad depth analysis results written to:", output_file, "\n")
  cat("Number of analyses:", nrow(result), "\n")
  cat("Coefficient used:", unique(result$coefficient), "\n")
  cat("Bootstrap iterations:", unique(result$n_bootstrap), "\n")
  cat("Confidence level:", unique(result$conf_level), "%\n")
  
  # Report significance counts
  if ("p_value_exceeds_threshold" %in% colnames(result)) {
    valid_p <- !is.na(result$p_value_exceeds_threshold)
    sig_raw <- sum(result$p_value_exceeds_threshold < 0.05, na.rm = TRUE)
    cat("Significant results (raw p < 0.05):", sig_raw, "/", sum(valid_p), "\n")
    
    if ("p_value_fdr" %in% colnames(result)) {
      sig_fdr <- sum(result$p_value_fdr < 0.05, na.rm = TRUE)
      cat("Significant results (FDR corrected p < 0.05):", sig_fdr, "/", sum(valid_p), "\n")
    }
  }
  
} else {
  result <- data.frame(
    types_of_interest = character(0),
    ref_superficial = character(0),
    ref_deep = character(0),
    n_types_interest = integer(0),
    n_ref_superficial = integer(0),
    n_ref_deep = integer(0),
    n_samples_interest = integer(0),
    n_samples_superficial = integer(0),
    n_samples_deep = integer(0),
    n_neurons_interest = integer(0),
    n_neurons_superficial = integer(0),
    n_neurons_deep = integer(0),
    observed_sup_d = numeric(0),
    observed_deep_d = numeric(0),
    observed_distance_diff = numeric(0),
    observed_delta_thres_base = numeric(0),
    observed_delta_thres = numeric(0),
    observed_test_statistic = numeric(0),
    observed_bias_ratio = numeric(0),
    bootstrap_distance_diff_median = numeric(0),
    bootstrap_distance_diff_lower = numeric(0),
    bootstrap_distance_diff_upper = numeric(0),
    bootstrap_delta_thres_base_median = numeric(0),
    bootstrap_delta_thres_base_lower = numeric(0),
    bootstrap_delta_thres_base_upper = numeric(0),
    bootstrap_test_stat_median = numeric(0),
    bootstrap_test_stat_lower = numeric(0),
    bootstrap_test_stat_upper = numeric(0),
    bootstrap_bias_ratio_median = numeric(0),
    bootstrap_bias_ratio_lower = numeric(0),
    bootstrap_bias_ratio_upper = numeric(0),
    coefficient = numeric(0),
    p_value_exceeds_threshold = numeric(0),
    p_value_bias_ratio = numeric(0),
    direction = character(0),
    conf_level = numeric(0),
    n_bootstrap = integer(0),
    p_value_fdr = numeric(0),
    significant_fdr = logical(0),
    split = character(0),
    neuropil = character(0),
    syn_type = character(0),
    preset = character(0)
  )
  cat("No valid analyses could be performed\n")
}
fwrite(result, output_file)
