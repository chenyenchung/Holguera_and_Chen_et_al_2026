#!/usr/bin/env Rscript
renv::load("/scratch/ycc520/flyem")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(Rcpp))

# Parse command-line arguments
argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

if (interactive()) {
  source("src/utils.r")
  argvs$np <- "LO_L"
  argvs$syn_type <- "pre"
  argvs$ts <- "data/selectors.csv"
  argvs$ann <- "data/visual_neurons_anno.csv"
  argvs$meta <- "data/viz_meta.csv"
  argvs$ref_groups <- "data/reference_groups.csv"
  argvs$sparse_limit <- 100L
  argvs$coefficient <- 0.5
  argvs$n_bootstrap <- 1000L
  argvs$conf_int <- 95.0
  argvs$cppsrc <- "src/stats/bin/broad_depth.cpp"
  argvs$batch_id <- "000"  # No batching in interactive mode
  argvs$genes <- NULL      # NULL = process all genes
  syn_path <- file.path("int/idv_mat/", paste0(argvs$np, "_rotated.csv.gz"))
} else {
  argvs$sparse_limit <- as.integer(argvs$sparse_limit)
  argvs$coefficient <- if (is.null(argvs$coefficient)) 0.5 else as.numeric(argvs$coefficient)
  argvs$n_bootstrap <- if (is.null(argvs$n_bootstrap)) 1000L else as.integer(argvs$n_bootstrap)
  argvs$conf_int <- if (is.null(argvs$conf_int)) 95.0 else as.numeric(argvs$conf_int)
  syn_path <- argvs$synf

  # Robust file path resolution for C++ source
  if (is.null(argvs$cppsrc)) {
    candidates <- c(
      "./broad_depth.cpp",
      "src/stats/bin/broad_depth.cpp"
    )
    found <- candidates[file.exists(candidates)][1]
    if (is.na(found)) {
      stop("Cannot find broad_depth.cpp in expected locations: ",
           paste(candidates, collapse = ", "))
    }
    argvs$cppsrc <- found
  }
}

# Verify C++ source file exists before sourcing
if (!file.exists(argvs$cppsrc)) {
  stop(sprintf("C++ source file not found: %s", argvs$cppsrc))
}

# Compile C++ source
cat("Compiling C++ bootstrap code...\n")
Rcpp::sourceCpp(argvs$cppsrc)

# Load and validate metadata
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

# Load annotations
opc_anno <- validate_config_file(argvs$ann, c("cell_type"))

# Load synapse coordinates
cat("Loading synapse coordinates...\n")
np_coord <- fread(syn_path)

# Store unfiltered copy for reference group calculation
np_coord_unfiltered <- copy(np_coord)

# Apply sparse filtering
cat("Filtering sparse cell types...\n")
np_coord <- filter_type(
  np_coord, syn_type = argvs$syn_type, sparse_limit = argvs$sparse_limit
)

# Load gene expression data and merge onto synapses
cat("Merging gene expression data...\n")
selectors <- fread(argvs$ts, header = TRUE)
np_coord <- filter_geneexp(np_coord, opc_anno, selectors, syn_type = argvs$syn_type)

print(head(np_coord))

# Get available cell types for reference group filtering
available_groups <- unique(
  c(np_coord_unfiltered$pre_type, np_coord_unfiltered$post_type)
)

# Load reference groups from configuration file
ref_groups_file <- argvs$ref_groups
if (!file.exists(ref_groups_file)) {
  stop(sprintf("Reference groups file not found: %s", ref_groups_file))
}

ref_groups_config <- fread(ref_groups_file)
ref_config <- ref_groups_config[neuropil == argvs$np]

if (nrow(ref_config) == 0) {
  stop(sprintf("No reference groups found for neuropil: %s", argvs$np))
}

# Parse base reference groups based on pattern type
if (ref_config$pattern_type == "regex") {
  ref_superficial_all <- grep(ref_config$ref_superficial, available_groups, value = TRUE)
  ref_deep_all <- grep(ref_config$ref_deep, available_groups, value = TRUE)
} else if (ref_config$pattern_type == "exact") {
  ref_superficial_all <- intersect(
    strsplit(ref_config$ref_superficial, ";")[[1]],
    available_groups
  )
  ref_deep_all <- intersect(
    strsplit(ref_config$ref_deep, ";")[[1]],
    available_groups
  )
} else {
  stop("Unknown pattern_type in reference_groups.csv: ", ref_config$pattern_type)
}

# # Create Notch-stratified reference groups
# # Get Notch annotations for reference types
# ref_notch_anno <- opc_anno[
#   cell_type %in% c(ref_superficial_all, ref_deep_all),
#   .(cell_type, Notch)
# ]
# 
# # Split reference groups by Notch category
# ref_superficial_all <- intersect(
#   ref_superficial_all,
#   ref_notch_anno[Notch == "Notch On", cell_type]
# )
# ref_superficial_all <- intersect(
#   ref_superficial_all,
#   ref_notch_anno[Notch == "Notch Off", cell_type]
# )
# 
# ref_deep_all <- intersect(
#   ref_deep_all,
#   ref_notch_anno[Notch == "Notch On", cell_type]
# )
# ref_deep_all <- intersect(
#   ref_deep_all,
#   ref_notch_anno[Notch == "Notch Off", cell_type]
# )

# Log reference group composition
cat(sprintf("\n=== Reference Groups for %s (syn_type=%s) ===\n", argvs$np, argvs$syn_type))
cat(sprintf("Shared superficial (%d types): %s\n",
            length(ref_superficial_all),
            paste(sort(ref_superficial_all), collapse = ", ")))
cat(sprintf("Shared deep (%d types): %s\n",
            length(ref_deep_all),
            paste(sort(ref_deep_all), collapse = ", ")))
cat("==========================================\n\n")

# Helper function to prepare neuron-to-synapse mapping for neuron-level bootstrap
# Adapted from broad_depth.r
prepare_neuron_depth_mapping <- function(data, neuron_types, syn_type, zinvert) {
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

  # Apply z-axis inversion if needed
  if (zinvert) {
    all_data$depth <- all_data$depth * -1
  }

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

# Prepare shared reference group neuron mappings
cat("Preparing shared reference group neuron mappings...\n")

data_superficial <- np_coord_unfiltered[
  pre_type %in% ref_superficial_all | post_type %in% ref_superficial_all
]
data_deep <- np_coord_unfiltered[
  pre_type %in% ref_deep_all | post_type %in% ref_deep_all
]

ref_sup_mapping <- prepare_neuron_depth_mapping(
  data_superficial, ref_superficial_all, argvs$syn_type, plot_meta$zinvert
)
ref_deep_mapping <- prepare_neuron_depth_mapping(
  data_deep, ref_deep_all, argvs$syn_type, plot_meta$zinvert
)

# Helper function to extract gene columns (logical, excluding annotation columns)
get_gene_columns <- function(data) {
  logical_cols <- names(data)[sapply(data, is.logical)]
  exclude <- c(
    "Confident_annotation", "newly_ann", "putative_OPC",
    grep("^putative_hl[0-9]+$", names(data), value = TRUE)
  )
  return(setdiff(logical_cols, exclude))
}

# Helper function to create NA result for skipped genes
create_na_result <- function(gene, reason, notch_category = NA_character_) {
  data.frame(
    types_of_interest = gene,
    notch_category = notch_category,
    skip_reason = reason,
    ref_superficial = NA_character_,
    ref_deep = NA_character_,
    n_neurons_expressing = NA_integer_,
    n_synapses_expressing = NA_integer_,
    n_types_expressing = NA_integer_,
    types_expressing = NA_character_,
    conflict_with_references = NA,
    overlap_superficial = NA_integer_,
    overlap_deep = NA_integer_,
    overlap_any = NA,
    observed_sup_d = NA_real_,
    observed_deep_d = NA_real_,
    observed_distance_diff = NA_real_,
    observed_delta_thres_base = NA_real_,
    observed_delta_thres = NA_real_,
    observed_test_statistic = NA_real_,
    observed_bias_ratio = NA_real_,
    bootstrap_distance_diff_median = NA_real_,
    bootstrap_distance_diff_lower = NA_real_,
    bootstrap_distance_diff_upper = NA_real_,
    bootstrap_delta_thres_base_median = NA_real_,
    bootstrap_delta_thres_base_lower = NA_real_,
    bootstrap_delta_thres_base_upper = NA_real_,
    bootstrap_test_stat_median = NA_real_,
    bootstrap_test_stat_lower = NA_real_,
    bootstrap_test_stat_upper = NA_real_,
    bootstrap_bias_ratio_median = NA_real_,
    bootstrap_bias_ratio_lower = NA_real_,
    bootstrap_bias_ratio_upper = NA_real_,
    p_value_exceeds_threshold = NA_real_,
    p_value_bias_ratio = NA_real_,
    direction = NA_character_,
    coefficient = argvs$coefficient,
    n_bootstrap = argvs$n_bootstrap,
    conf_level = argvs$conf_int,
    stringsAsFactors = FALSE
  )
}

# Extract gene columns
all_gene_symbols <- get_gene_columns(np_coord)

# Filter to batch subset if specified
if (!is.null(argvs$genes) && argvs$genes != "") {
  requested_genes <- strsplit(argvs$genes, ",")[[1]]
  # Validate requested genes exist
  missing_genes <- setdiff(requested_genes, all_gene_symbols)
  if (length(missing_genes) > 0) {
    warning(sprintf("Requested genes not found: %s", paste(missing_genes, collapse = ", ")))
  }
  gene_symbols <- intersect(requested_genes, all_gene_symbols)
  cat(sprintf("Processing batch %s: %d genes (of %d total)\n",
              argvs$batch_id, length(gene_symbols), length(all_gene_symbols)))
} else {
  gene_symbols <- all_gene_symbols
  cat(sprintf("Processing all genes: %d\n", length(gene_symbols)))
}

if (length(gene_symbols) == 0) {
  stop("No valid genes to process in this batch")
}

# Gene loop: test each gene separately
results <- list()
n_tested <- 0
n_skipped <- 0

cat("\n=== Starting Gene-Level Bootstrap Analysis (Notch-Stratified) ===\n")

for (gene in gene_symbols) {
  cat(sprintf("Processing gene: %s...\n", gene))

  # Define Notch categories to test
  notch_categories <- c("Notch On", "Notch Off")

  for (notch_cat in notch_categories) {
    cat(sprintf("  Notch category: %s... ", notch_cat))

    # STEP 1: Filter to expressing synapses IN THIS NOTCH CATEGORY
    expressing_synapses <- np_coord[get(gene) == TRUE & Notch == notch_cat]

    # STEP 2: Edge case checks

    # Edge case: no expressing synapses in this Notch category
    if (nrow(expressing_synapses) == 0) {
      cat("no expressing neurons in this category\n")
      result_key <- paste(gene, notch_cat, sep = "__")
      results[[result_key]] <- create_na_result(
        gene, "no_expressing_neurons", notch_category = notch_cat
      )
      n_skipped <- n_skipped + 1
      next
    }

    # Edge case: insufficient synapses
    if (nrow(expressing_synapses) < argvs$sparse_limit) {
      cat(sprintf("insufficient synapses (%d < %d)\n",
                  nrow(expressing_synapses), argvs$sparse_limit))
      result_key <- paste(gene, notch_cat, sep = "__")
      results[[result_key]] <- create_na_result(
        gene, "insufficient_synapses", notch_category = notch_cat
      )
      n_skipped <- n_skipped + 1
      next
    }

    # Check neuron count
    neuron_col <- paste0(argvs$syn_type, "_root_id")
    n_neurons <- expressing_synapses[, uniqueN(get(neuron_col))]

    # Edge case: insufficient neurons
    if (n_neurons < 3) {
      cat(sprintf("insufficient neurons (%d < 3)\n", n_neurons))
      result_key <- paste(gene, notch_cat, sep = "__")
      results[[result_key]] <- create_na_result(
        gene, "insufficient_neurons", notch_category = notch_cat
      )
      n_skipped <- n_skipped + 1
      next
    }

    # Get expressing cell types
    type_col <- paste0(argvs$syn_type, "_type")
    expressing_types <- unique(expressing_synapses[[type_col]])

    # STEP 3: Use shared reference groups/mappings for all Notch strata
    ref_sup_current <- ref_superficial_all
    ref_deep_current <- ref_deep_all

    # Check overlap with reference groups (references are shared across Notch strata)
    overlap_sup <- sum(expressing_types %in% ref_sup_current)
    overlap_deep <- sum(expressing_types %in% ref_deep_current)
    overlap_any <- (overlap_sup > 0 || overlap_deep > 0)
    conflict_flag <- (overlap_sup > 0 && overlap_deep > 0)

    # STEP 4: Prepare interest group mapping
    interest_data <- expressing_synapses[, .(
      neuron_id = get(neuron_col),
      depth = get(paste0(argvs$syn_type, "_rz"))
    )]

    # Apply z-axis inversion if needed
    if (plot_meta$zinvert) {
      interest_data$depth <- interest_data$depth * -1
    }

    setorder(interest_data, neuron_id)
    interest_neuron_info <- interest_data[, .N, by = neuron_id]
    interest_neuron_info[, start_idx := c(0, cumsum(N[-.N]))]

    interest_mapping <- list(
      depths = interest_data$depth,
      neuron_starts = interest_neuron_info$start_idx,
      neuron_counts = interest_neuron_info$N
    )

    # STEP 5: Run bootstrap (C++ function - unchanged)
    bootstrap_result <- tryCatch(
      {
        perform_broad_depth_bootstrap_neuron_level(
          ref_sup_depths = ref_sup_mapping$depths,
          ref_sup_neuron_starts = as.integer(ref_sup_mapping$neuron_starts),
          ref_sup_neuron_counts = as.integer(ref_sup_mapping$neuron_counts),
          ref_deep_depths = ref_deep_mapping$depths,
          ref_deep_neuron_starts = as.integer(ref_deep_mapping$neuron_starts),
          ref_deep_neuron_counts = as.integer(ref_deep_mapping$neuron_counts),
          interest_depths = interest_mapping$depths,
          interest_neuron_starts = as.integer(interest_mapping$neuron_starts),
          interest_neuron_counts = as.integer(interest_mapping$neuron_counts),
          coefficient = argvs$coefficient,
          n_bootstrap = argvs$n_bootstrap,
          conf_int = argvs$conf_int,
          seed = NULL
        )
      },
      error = function(e) {
        cat(sprintf("ERROR: %s\n", e$message))
        return(NULL)
      }
    )

    # Handle bootstrap failure
    if (is.null(bootstrap_result)) {
      result_key <- paste(gene, notch_cat, sep = "__")
      results[[result_key]] <- create_na_result(
        gene, "bootstrap_failed", notch_category = notch_cat
      )
      n_skipped <- n_skipped + 1
      next
    }

    # Determine direction
    if (abs(bootstrap_result$observed_distance_diff) > bootstrap_result$observed_delta_thres) {
      direction <- if (bootstrap_result$observed_distance_diff > 0) "superficial" else "deep"
    } else {
      direction <- "neither"
    }

    cat(sprintf("OK (%d neurons, %d synapses, direction=%s)\n",
                n_neurons, nrow(expressing_synapses), direction))

    # STEP 6: Store result with notch_category column
    result_key <- paste(gene, notch_cat, sep = "__")
    results[[result_key]] <- data.frame(
      types_of_interest = gene,
      notch_category = notch_cat,
      skip_reason = NA_character_,
      ref_superficial = paste(ref_sup_current, collapse = ";"),
      ref_deep = paste(ref_deep_current, collapse = ";"),
    n_neurons_expressing = n_neurons,
    n_synapses_expressing = nrow(expressing_synapses),
    n_types_expressing = length(expressing_types),
    types_expressing = paste(sort(expressing_types), collapse = ";"),
    conflict_with_references = conflict_flag,
    overlap_superficial = overlap_sup,
    overlap_deep = overlap_deep,
    overlap_any = overlap_any,
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
    p_value_exceeds_threshold = bootstrap_result$p_value_exceeds_threshold,
    p_value_bias_ratio = bootstrap_result$p_value_bias_ratio,
    direction = direction,
    coefficient = argvs$coefficient,
    n_bootstrap = argvs$n_bootstrap,
    conf_level = argvs$conf_int,
      stringsAsFactors = FALSE
    )

    n_tested <- n_tested + 1
  }  # End Notch category loop
}  # End gene loop

# Combine all results
result_df <- rbindlist(results, fill = TRUE)

# Apply provisional within-job FDR correction stratified by Notch category.
# combine_selector_results.r replaces these values after all workflow batches
# are collected so the final result is independent of genes_per_batch.
result_df$p_value_fdr <- NA_real_
result_df$significant_fdr <- NA

for (notch_cat in c("Notch On", "Notch Off")) {
  subset_idx <- result_df$notch_category == notch_cat
  valid_p <- subset_idx & !is.na(result_df$p_value_exceeds_threshold)

  if (sum(valid_p) > 0) {
    result_df$p_value_fdr[valid_p] <- p.adjust(
      result_df$p_value_exceeds_threshold[valid_p],
      method = "fdr"
    )
    result_df$significant_fdr[valid_p] <- result_df$p_value_fdr[valid_p] < 0.05
  }
}

# Add metadata
result_df$neuropil <- argvs$np
result_df$syn_type <- argvs$syn_type
result_df$batch_id <- if (is.null(argvs$batch_id)) "000" else argvs$batch_id
result_df$sparse_limit <- argvs$sparse_limit
result_df$analysis_date <- as.character(Sys.Date())

# Write output
if (is.null(argvs$batch_id) || argvs$batch_id == "000") {
  # No batching - backward compatible filename
  output_file <- paste0(argvs$np, "_", argvs$syn_type, "_selector_depth.csv")
} else {
  # Batching enabled - include batch ID
  output_file <- paste0(argvs$np, "_", argvs$syn_type, "_batch", argvs$batch_id, "_selector_depth.csv")
}
fwrite(result_df, output_file)

# Print summary
cat("\n=== Selector Depth Analysis Summary ===\n")
cat("Neuropil:", argvs$np, "\n")
cat("Syn type:", argvs$syn_type, "\n")
cat("Total analyses:", nrow(result_df), "\n")
cat("Unique genes:", length(unique(result_df$types_of_interest)), "\n")
cat("Genes tested:", n_tested, "\n")
cat("Genes skipped:", n_skipped, "\n\n")

for (notch_cat in c("Notch On", "Notch Off")) {
  subset <- result_df[notch_category == notch_cat]
  valid_p <- !is.na(subset$p_value_exceeds_threshold)

  cat(sprintf("=== %s ===\n", notch_cat))
  cat("Total tests:", nrow(subset), "\n")

  if (sum(valid_p) > 0) {
    n_sig_raw <- sum(subset$p_value_exceeds_threshold < 0.05, na.rm = TRUE)
    n_sig_fdr <- sum(subset$significant_fdr, na.rm = TRUE)
    n_superficial <- sum(subset$direction == "superficial", na.rm = TRUE)
    n_deep <- sum(subset$direction == "deep", na.rm = TRUE)
    n_neither <- sum(subset$direction == "neither", na.rm = TRUE)

    cat("Significant (raw p < 0.05):", n_sig_raw, "/", sum(valid_p), "\n")
    cat("Significant (provisional within-job FDR < 0.05):",
        n_sig_fdr, "/", sum(valid_p), "\n")
    cat("Direction: superficial=", n_superficial,
        ", deep=", n_deep,
        ", neither=", n_neither, "\n\n")
  } else {
    cat("No valid tests in this category\n\n")
  }
}

if (sum(result_df$conflict_with_references, na.rm = TRUE) > 0) {
  cat("WARNING: Genes with reference conflicts:",
      sum(result_df$conflict_with_references, na.rm = TRUE), "\n")
}

cat("\nOutput written to:", output_file, "\n")
cat("========================================\n")
