#!/usr/bin/env Rscript
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

required_args <- c("np", "synf", "ann", "meta", "ref_groups")
missing_args <- required_args[vapply(required_args, function(x) is.null(argvs[[x]]), logical(1))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

argvs$sparse_limit <- if (is.null(argvs$sparse_limit)) 100L else as.integer(argvs$sparse_limit)
argvs$min_neurons <- if (is.null(argvs$min_neurons)) 3L else as.integer(argvs$min_neurons)
argvs$type_chunk_size <- if (is.null(argvs$type_chunk_size)) 5L else as.integer(argvs$type_chunk_size)
if (argvs$type_chunk_size < 1L) stop("--type_chunk_size must be >= 1")
argvs$syn_types <- if (is.null(argvs$syn_types)) "pre,post" else argvs$syn_types
syn_types <- trimws(strsplit(argvs$syn_types, "[,;]")[[1]])
syn_types <- syn_types[nzchar(syn_types)]
unknown_syn_types <- setdiff(syn_types, c("pre", "post"))
if (length(unknown_syn_types) > 0) {
  stop("Unknown syn_types: ", paste(unknown_syn_types, collapse = ", "))
}

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

interest <- list()
for (syn_type in syn_types) {
  type_col <- paste0(syn_type, "_type")
  root_col <- paste0(syn_type, "_root_id")
  depth_col <- paste0(syn_type, "_rz")
  interest[[syn_type]] <- list()

  for (cell_type in eligible$cell_type) {
    type_syn <- syn[get(type_col) == cell_type]
    n_syn <- nrow(type_syn)
    n_neurons <- if (n_syn > 0) type_syn[, uniqueN(get(root_col))] else 0L
    n_depths <- if (n_syn > 0) type_syn[, uniqueN(get(depth_col), na.rm = TRUE)] else 0L

    skip_reason <- NA_character_
    mapping <- NULL
    if (n_syn < argvs$sparse_limit) {
      skip_reason <- "insufficient_synapses"
    } else if (n_neurons < argvs$min_neurons) {
      skip_reason <- "insufficient_neurons"
    } else if (n_depths < 2) {
      skip_reason <- "insufficient_depth_variation"
    } else {
      mapping <- prepare_interest_mapping(type_syn, syn_type, zinvert)
      if (length(mapping$depths) == 0) {
        skip_reason <- "no_depths"
        mapping <- NULL
      }
    }

    interest[[syn_type]][[cell_type]] <- list(
      skip_reason = skip_reason,
      mapping = mapping
    )
  }
}

prepared <- list(
  np = argvs$np,
  eligible = eligible,
  refs = refs,
  ref_sup_mapping = ref_sup_mapping,
  ref_deep_mapping = ref_deep_mapping,
  interest = interest
)

prepared_out <- sprintf("%s_type_depth_prepared.rds", argvs$np)
chunks_out <- sprintf("%s_type_depth_chunks.csv", argvs$np)
saveRDS(prepared, prepared_out)

chunks <- rbindlist(lapply(syn_types, function(syn_type) {
  cell_types <- eligible$cell_type
  shard_ids <- ceiling(seq_along(cell_types) / argvs$type_chunk_size)
  data.table(
    np = argvs$np,
    syn_type = syn_type,
    shard_id = shard_ids,
    cell_types = cell_types
  )[, .(cell_types = paste(cell_types, collapse = ";")), by = .(np, syn_type, shard_id)]
}))
fwrite(chunks, chunks_out)

cat(sprintf(
  "Wrote %s and %s with %d shard rows\n",
  prepared_out, chunks_out, nrow(chunks)
))
