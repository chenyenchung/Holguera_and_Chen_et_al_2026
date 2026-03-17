#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

### TODO: Interactive mode for development
if (interactive()) {
  source("./src/utils.r", chdir = FALSE)
  argvs$output <- "partner_summary.csv"
  # Paths to idv_mat files
  argvs$me_l <- "int/idv_mat/ME_L_rotated.csv.gz"
  argvs$me_r <- "int/idv_mat/ME_R_rotated.csv.gz"
  argvs$lo_l <- "int/idv_mat/LO_L_rotated.csv.gz"
  argvs$lo_r <- "int/idv_mat/LO_R_rotated.csv.gz"
  argvs$lop_l <- "int/idv_mat/LOP_L_rotated.csv.gz"
  argvs$lop_r <- "int/idv_mat/LOP_R_rotated.csv.gz"
}

# Function to process a single neuropil file
process_neuropil_file <- function(file_path, neuropil_name) {
  cat(sprintf("Processing %s...\n", basename(file_path)))
  data <- fread(file_path, colClasses = "character")
  
  out <- data[
    , .N, by = c("pre_root_id", "post_root_id", "pre_type", "post_type")
  ][pre_type != "" & post_type != ""]
  
  # Consider > 4 synapses per pair to be real
  to_keep <- out[
   , .(nSynapse = sum(N)),  by = c("pre_root_id", "post_root_id")
  ][nSynapse > 4]
  
  out <- merge(out, to_keep, by = c("pre_root_id", "post_root_id"))
  out <- out[
    , .(neuropil = neuropil_name, nSynapse = sum(N), nNeuron = .N),
    by = c("pre_type", "post_type")
  ]
  return(out)
}

# Function to combine bilateral neuropil data
combine_bilateral_neuropil <- function(left_file, right_file, neuropil_name) {
  cat(sprintf("Combining bilateral data for %s...\n", neuropil_name))
  
  # Process each file
  left_data <- process_neuropil_file(left_file, neuropil_name)
  right_data <- process_neuropil_file(right_file, neuropil_name)
  
  out <- rbind(left_data, right_data)
  out <- out[
    , .(
      neuropil = neuropil_name,
      nSynapse = sum(nSynapse), nNeuron = sum(nNeuron)
    ),
    by = c("pre_type", "post_type")
  ]
  
  return(out)
}

cat("Starting partner extraction...\n")

# Process each neuropil
neuropils <- list(
  "ME" = list(left = argvs$me_l, right = argvs$me_r),
  "LO" = list(left = argvs$lo_l, right = argvs$lo_r),
  "LOP" = list(left = argvs$lop_l, right = argvs$lop_r)
)

all_data <- list()

for (neuropil_name in names(neuropils)) {
  neuropil_data <- combine_bilateral_neuropil(
    neuropils[[neuropil_name]]$left,
    neuropils[[neuropil_name]]$right,
    neuropil_name
  )
  all_data[[neuropil_name]] <- neuropil_data
}

# Combine all neuropil data
final_data <- do.call(rbind, all_data)

# Reorder columns to match specification: neuropil, pre, post, nSynapse, nNeuron
final_data <- final_data[ , .(neuropil, pre_type, post_type, nSynapse, nNeuron)]
# Write output
fwrite(final_data, argvs$output, row.names = FALSE)

cat(sprintf("Partner extraction complete! Output saved to: %s\n", argvs$output))
cat(sprintf("Total rows: %d\n", nrow(final_data)))
cat(sprintf("Neuropils: %s\n", paste(unique(final_data$neuropil), collapse = ", ")))
