#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

if (interactive()) {
  argvs$syn <- "work/6c/fd9678acb41db6ff64e7b4e93cb2c9/ME_L.csv.gz"
}

if (is.null(argvs$syn) || !file.exists(argvs$syn)) {
  stop("Cannot find the csv file containing synases.")
}

np <- fread(argvs$syn)
np_coord <- rbind(
  np[, .(x = pre_x, y = pre_y, z = pre_z)],
  np[, .(x = post_x, y = post_y, z = post_z)]
)

rotate <- prcomp(np_coord)$x

np[, pre_rx := rotate[seq_len(nrow(np)), "PC1"]]
np[, pre_ry := rotate[seq_len(nrow(np)), "PC2"]]
np[, pre_rz := rotate[seq_len(nrow(np)), "PC3"]]

next_half <- seq(nrow(np) + 1, 2 * nrow(np))
np[, post_rx := rotate[next_half, "PC1"]]
np[, post_ry := rotate[next_half, "PC2"]]
np[, post_rz := rotate[next_half, "PC3"]]

outname <- sub("\\.csv\\.gz$", "", argvs$syn)

fwrite(np, paste0(outname, "_rotated.csv.gz"))
