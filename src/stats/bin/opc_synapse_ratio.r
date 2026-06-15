#!/usr/bin/env Rscript
renv::load("/scratch/ycc520/flyem")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

require_arg <- function(name) {
  value <- argvs[[name]]
  if (is.null(value) || is.na(value) || value == "") {
    stop(sprintf("Missing required argument: --%s", name))
  }
  value
}

validate_columns <- function(data, required_cols, label) {
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "%s is missing required columns: %s",
      label,
      paste(missing_cols, collapse = ", ")
    ))
  }
}

is_true <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

analyze_synapse_file <- function(syn_path, neuropil, opc_types) {
  syn <- fread(syn_path, select = c("pre_type", "post_type"))
  validate_columns(syn, c("pre_type", "post_type"), syn_path)

  total_synapses <- nrow(syn)
  side <- sub("^.*_", "", neuropil)

  data.table(
    neuropil = neuropil,
    side = side,
    syn_type = c("pre", "post"),
    opc_synapse_count = c(
      sum(syn$pre_type %chin% opc_types),
      sum(syn$post_type %chin% opc_types)
    ),
    total_synapse_count = total_synapses,
    n_putative_opc_types = length(opc_types)
  )[
    ,
    opc_synapse_ratio := fifelse(
      total_synapse_count > 0,
      opc_synapse_count / total_synapse_count,
      NA_real_
    )
  ][
    ,
    .(
      neuropil,
      side,
      syn_type,
      opc_synapse_count,
      total_synapse_count,
      opc_synapse_ratio,
      n_putative_opc_types
    )
  ]
}

ann_path <- require_arg("ann")
lo_l_path <- require_arg("lo_l")
lo_r_path <- require_arg("lo_r")
output_file <- if (is.null(argvs$output) || argvs$output == "") {
  "LO_opc_synapse_ratio.csv"
} else {
  argvs$output
}

anno <- fread(ann_path)
validate_columns(anno, c("cell_type", "putative_OPC"), ann_path)

opc_types <- anno[is_true(putative_OPC), unique(cell_type)]
if (length(opc_types) == 0) {
  stop(sprintf("No putative OPC cell types found in %s", ann_path))
}

results <- rbindlist(list(
  analyze_synapse_file(lo_l_path, "LO_L", opc_types),
  analyze_synapse_file(lo_r_path, "LO_R", opc_types)
))

fwrite(results, output_file)

cat("OPC synapse ratio results written to:", output_file, "\n")
cat("Putative OPC types:", length(opc_types), "\n")
cat("Rows:", nrow(results), "\n")
