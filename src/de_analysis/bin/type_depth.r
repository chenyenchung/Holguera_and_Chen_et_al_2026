#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rcpp))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

required_args <- c("np", "synf", "syn_type", "ann", "meta", "ref_groups", "cppsrc")
missing_args <- required_args[vapply(required_args, function(x) is.null(argvs[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

argvs$sparse_limit <- if (is.null(argvs$sparse_limit)) 100L else as.integer(argvs$sparse_limit)
argvs$min_neurons <- if (is.null(argvs$min_neurons)) 3L else as.integer(argvs$min_neurons)
argvs$coefficient <- if (is.null(argvs$coefficient)) 0.5 else as.numeric(argvs$coefficient)
argvs$n_bootstrap <- if (is.null(argvs$n_bootstrap)) 1000L else as.integer(argvs$n_bootstrap)
argvs$conf_int <- if (is.null(argvs$conf_int)) 95 else as.numeric(argvs$conf_int)

validate_config_file <- function(file_path, required_cols) {
  if (!file.exists(file_path)) stop(sprintf("File not found: %s", file_path))
  data <- fread(file_path)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in %s: %s", file_path, paste(missing_cols, collapse = ", ")))
  }
  data
}

parse_reference_types <- function(ref_config, available_groups) {
  if (nrow(ref_config) != 1) {
    stop(sprintf("Expected exactly one reference row, found %d", nrow(ref_config)))
  }

  if (ref_config$pattern_type == "regex") {
    ref_superficial <- grep(ref_config$ref_superficial, available_groups, value = TRUE)
    ref_deep <- grep(ref_config$ref_deep, available_groups, value = TRUE)
  } else if (ref_config$pattern_type == "exact") {
    ref_superficial <- intersect(
      trimws(strsplit(ref_config$ref_superficial, ";")[[1]]),
      available_groups
    )
    ref_deep <- intersect(
      trimws(strsplit(ref_config$ref_deep, ";")[[1]]),
      available_groups
    )
  } else {
    stop("Unknown reference pattern_type: ", ref_config$pattern_type)
  }

  overlap <- intersect(ref_superficial, ref_deep)
  if (length(overlap) > 0) {
    stop(sprintf("Reference groups overlap for %s: %s", ref_config$neuropil, paste(overlap, collapse = ", ")))
  }

  list(superficial = sort(ref_superficial), deep = sort(ref_deep))
}

prepare_neuron_depth_mapping <- function(data, neuron_types, zinvert) {
  pre_part <- data[pre_type %in% neuron_types, .(
    neuron_id = paste0("pre_", pre_root_id),
    depth = pre_rz
  )]
  post_part <- data[post_type %in% neuron_types, .(
    neuron_id = paste0("post_", post_root_id),
    depth = post_rz
  )]
  all_data <- rbind(pre_part, post_part)
  all_data <- all_data[!is.na(depth)]
  if (zinvert) all_data[, depth := depth * -1]

  setorder(all_data, neuron_id)
  neuron_info <- all_data[, .N, by = neuron_id]
  neuron_info[, start_idx := c(0L, cumsum(N[-.N]))]

  list(
    depths = all_data$depth,
    neuron_starts = neuron_info$start_idx,
    neuron_counts = neuron_info$N
  )
}

prepare_interest_mapping <- function(data, syn_type, zinvert) {
  neuron_col <- paste0(syn_type, "_root_id")
  depth_col <- paste0(syn_type, "_rz")
  interest <- data[, .(
    neuron_id = get(neuron_col),
    depth = get(depth_col)
  )]
  interest <- interest[!is.na(depth)]
  if (zinvert) interest[, depth := depth * -1]

  setorder(interest, neuron_id)
  neuron_info <- interest[, .N, by = neuron_id]
  neuron_info[, start_idx := c(0L, cumsum(N[-.N]))]

  list(
    depths = interest$depth,
    neuron_starts = neuron_info$start_idx,
    neuron_counts = neuron_info$N
  )
}

make_result <- function(row, skip_reason = NA_character_, values = list()) {
  base <- data.table(
    cell_type = row$cell_type,
    ozel2021_cluster = row$ozel2021_cluster,
    Notch = row$Notch,
    ntype = row$ntype,
    func = row$func,
    temporal_label = row$temporal_label,
    neuropil = argvs$np,
    syn_type = argvs$syn_type,
    skip_reason = skip_reason,
    direction = NA_character_,
    significant_fdr = NA,
    p_value_fdr = NA_real_,
    p_value_exceeds_threshold = NA_real_,
    p_value_bias_ratio = NA_real_,
    observed_bias_ratio = NA_real_,
    observed_sup_d = NA_real_,
    observed_deep_d = NA_real_,
    observed_distance_diff = NA_real_,
    observed_delta_thres_base = NA_real_,
    observed_delta_thres = NA_real_,
    observed_test_statistic = NA_real_,
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
    n_neurons_interest = NA_integer_,
    n_syn_interest = NA_integer_,
    n_neurons_superficial = NA_integer_,
    n_syn_superficial = NA_integer_,
    n_neurons_deep = NA_integer_,
    n_syn_deep = NA_integer_,
    ref_superficial = paste(refs$superficial, collapse = ";"),
    ref_deep = paste(refs$deep, collapse = ";"),
    coefficient = argvs$coefficient,
    n_bootstrap = argvs$n_bootstrap,
    conf_level = argvs$conf_int,
    sparse_limit = argvs$sparse_limit,
    min_neurons = argvs$min_neurons
  )

  for (nm in names(values)) base[, (nm) := values[[nm]]]
  base
}

if (!file.exists(argvs$cppsrc)) stop("C++ source file not found: ", argvs$cppsrc)
Rcpp::sourceCpp(argvs$cppsrc)

plot_meta <- validate_config_file(
  argvs$meta,
  c("neuropil", "zinvert")
)[neuropil == argvs$np]
if (nrow(plot_meta) != 1) stop("Expected one metadata row for neuropil: ", argvs$np)
zinvert <- as.logical(plot_meta$zinvert[[1]])

ann <- validate_config_file(
  argvs$ann,
  c("cell_type", "Confident_annotation", "ozel2021_cluster", "Notch", "ntype", "func", "temporal_label")
)
eligible <- ann[Confident_annotation == "Y" & !is.na(ozel2021_cluster)]
eligible[, ozel2021_cluster := as.character(ozel2021_cluster)]
setorder(eligible, cell_type)

syn <- fread(argvs$synf)
required_syn_cols <- c(
  "pre_type", "post_type", "pre_root_id", "post_root_id",
  "pre_rz", "post_rz"
)
missing_syn_cols <- setdiff(required_syn_cols, names(syn))
if (length(missing_syn_cols) > 0) {
  stop("Missing synapse columns: ", paste(missing_syn_cols, collapse = ", "))
}

ref_config <- validate_config_file(
  argvs$ref_groups,
  c("neuropil", "ref_superficial", "ref_deep", "pattern_type")
)[neuropil == argvs$np]
available_groups <- unique(c(syn$pre_type, syn$post_type))
refs <- parse_reference_types(ref_config, available_groups)
if (length(refs$superficial) == 0 || length(refs$deep) == 0) {
  stop("Reference groups are empty for ", argvs$np)
}

ref_sup_data <- syn[pre_type %in% refs$superficial | post_type %in% refs$superficial]
ref_deep_data <- syn[pre_type %in% refs$deep | post_type %in% refs$deep]
ref_sup_mapping <- prepare_neuron_depth_mapping(ref_sup_data, refs$superficial, zinvert)
ref_deep_mapping <- prepare_neuron_depth_mapping(ref_deep_data, refs$deep, zinvert)

if (length(ref_sup_mapping$depths) == 0 || length(ref_deep_mapping$depths) == 0) {
  stop("Reference mappings are empty for ", argvs$np)
}

type_col <- paste0(argvs$syn_type, "_type")
root_col <- paste0(argvs$syn_type, "_root_id")
depth_col <- paste0(argvs$syn_type, "_rz")

results <- vector("list", nrow(eligible))
for (idx in seq_len(nrow(eligible))) {
  row <- eligible[idx]
  type_syn <- syn[get(type_col) == row$cell_type]

  if (nrow(type_syn) < argvs$sparse_limit) {
    results[[idx]] <- make_result(row, "insufficient_synapses")
    next
  }

  n_neurons <- type_syn[, uniqueN(get(root_col))]
  if (n_neurons < argvs$min_neurons) {
    results[[idx]] <- make_result(row, "insufficient_neurons")
    next
  }

  if (type_syn[, uniqueN(get(depth_col), na.rm = TRUE)] < 2) {
    results[[idx]] <- make_result(row, "insufficient_depth_variation")
    next
  }

  interest_mapping <- prepare_interest_mapping(type_syn, argvs$syn_type, zinvert)
  if (length(interest_mapping$depths) == 0) {
    results[[idx]] <- make_result(row, "no_depths")
    next
  }

  bootstrap_result <- tryCatch(
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
    ),
    error = function(e) e
  )

  if (inherits(bootstrap_result, "error")) {
    results[[idx]] <- make_result(row, "bootstrap_failed")
    next
  }

  direction <- if (abs(bootstrap_result$observed_distance_diff) > bootstrap_result$observed_delta_thres) {
    if (bootstrap_result$observed_distance_diff > 0) "superficial" else "deep"
  } else {
    "neither"
  }

  values <- list(
    skip_reason = NA_character_,
    direction = direction,
    p_value_exceeds_threshold = bootstrap_result$p_value_exceeds_threshold,
    p_value_bias_ratio = bootstrap_result$p_value_bias_ratio,
    observed_bias_ratio = bootstrap_result$observed_bias_ratio,
    observed_sup_d = bootstrap_result$observed_sup_d,
    observed_deep_d = bootstrap_result$observed_deep_d,
    observed_distance_diff = bootstrap_result$observed_distance_diff,
    observed_delta_thres_base = bootstrap_result$observed_delta_thres_base,
    observed_delta_thres = bootstrap_result$observed_delta_thres,
    observed_test_statistic = bootstrap_result$observed_test_statistic,
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
    n_neurons_interest = bootstrap_result$n_neurons_interest,
    n_syn_interest = bootstrap_result$n_syn_interest,
    n_neurons_superficial = bootstrap_result$n_neurons_sup,
    n_syn_superficial = bootstrap_result$n_syn_sup,
    n_neurons_deep = bootstrap_result$n_neurons_deep,
    n_syn_deep = bootstrap_result$n_syn_deep
  )
  results[[idx]] <- make_result(row, values = values)
}

result <- rbindlist(results, fill = TRUE)
valid <- !is.na(result$p_value_exceeds_threshold)
if (sum(valid) > 0) {
  result[valid, p_value_fdr := p.adjust(p_value_exceeds_threshold, method = "fdr")]
  result[valid, significant_fdr := p_value_fdr < 0.05]
}

output_file <- sprintf("%s_%s_type_depth.csv", argvs$np, argvs$syn_type)
fwrite(result, output_file)

cat(sprintf("Wrote %s with %d rows (%d tested)\n", output_file, nrow(result), sum(valid)))
