#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValue = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

if (is.null(argvs$syn) || !file.exists(argvs$syn)) {
  stop("Cannot find the csv file containing synases.")
}

syn_coord <- fread(argvs$syn)
syn_coord <- split(syn_coord, syn_coord$neuropil)

for (npil in names(syn_coord)) {
  message("Saving ", npil, "...")
  fwrite(syn_coord[[npil]], paste0(npil, ".csv.gz"))
}
