#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Rcpp))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

required_args <- c("np", "prepared", "syn_type", "shard_id", "cell_types", "cppsrc")
missing_args <- required_args[vapply(required_args, function(x) is.null(argvs[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

argvs$shard_id <- as.integer(argvs$shard_id)
argvs$sparse_limit <- if (is.null(argvs$sparse_limit)) 100L else as.integer(argvs$sparse_limit)
argvs$min_neurons <- if (is.null(argvs$min_neurons)) 3L else as.integer(argvs$min_neurons)
argvs$coefficient <- if (is.null(argvs$coefficient)) 0.5 else as.numeric(argvs$coefficient)
argvs$n_bootstrap <- if (is.null(argvs$n_bootstrap)) 1000L else as.integer(argvs$n_bootstrap)
argvs$conf_int <- if (is.null(argvs$conf_int)) 95 else as.numeric(argvs$conf_int)
argvs$bootstrap_seed <- if (is.null(argvs$bootstrap_seed)) 1L else as.integer(argvs$bootstrap_seed)

cell_types <- trimws(strsplit(argvs$cell_types, ";", fixed = TRUE)[[1]])
cell_types <- cell_types[nzchar(cell_types)]
if (length(cell_types) == 0) stop("No cell types supplied for shard")

if (!file.exists(argvs$cppsrc)) stop("C++ source file not found: ", argvs$cppsrc)
Rcpp::sourceCpp(argvs$cppsrc)

prepared <- readRDS(argvs$prepared)
if (!identical(prepared$np, argvs$np)) {
  stop(sprintf("Prepared input neuropil mismatch: expected %s, found %s", argvs$np, prepared$np))
}
if (!argvs$syn_type %in% names(prepared$interest)) {
  stop("Unknown syn_type in prepared input: ", argvs$syn_type)
}

refs <- prepared$refs
ref_sup_mapping <- prepared$ref_sup_mapping
ref_deep_mapping <- prepared$ref_deep_mapping
eligible <- copy(prepared$eligible)
setkey(eligible, cell_type)

stable_seed <- function(np, syn_type, cell_type, base_seed) {
  raw <- charToRaw(paste(np, syn_type, cell_type, base_seed, sep = "|"))
  value <- as.numeric(base_seed)
  for (byte in as.integer(raw)) {
    value <- (value * 1103515245 + byte + 12345) %% 2147483647
  }
  as.integer(max(1, value))
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

run_one <- function(cell_type) {
  row <- eligible[cell_type]
  if (nrow(row) != 1) stop("Unknown or duplicated cell_type in shard: ", cell_type)

  prepared_interest <- prepared$interest[[argvs$syn_type]][[cell_type]]
  if (is.null(prepared_interest)) {
    return(make_result(row, "no_prepared_mapping"))
  }

  if (!is.na(prepared_interest$skip_reason)) {
    return(make_result(row, prepared_interest$skip_reason))
  }

  interest_mapping <- prepared_interest$mapping
  if (is.null(interest_mapping) || length(interest_mapping$depths) == 0) {
    return(make_result(row, "no_depths"))
  }

  seed <- stable_seed(argvs$np, argvs$syn_type, cell_type, argvs$bootstrap_seed)
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
      seed = seed
    ),
    error = function(e) e
  )

  if (inherits(bootstrap_result, "error")) {
    return(make_result(row, "bootstrap_failed"))
  }

  direction <- if (abs(bootstrap_result$observed_distance_diff) > bootstrap_result$observed_delta_thres) {
    if (bootstrap_result$observed_distance_diff > 0) "superficial" else "deep"
  } else {
    "neither"
  }

  make_result(row, values = list(
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
  ))
}

result <- rbindlist(lapply(cell_types, run_one), fill = TRUE)

output_file <- sprintf(
  "%s_%s_shard%03d_type_depth.csv",
  argvs$np,
  argvs$syn_type,
  argvs$shard_id
)
fwrite(result, output_file)

cat(sprintf("Wrote %s with %d rows\n", output_file, nrow(result)))
