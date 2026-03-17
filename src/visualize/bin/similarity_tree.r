#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Cairo))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(tidytree))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(R.utils))

# Parse command line arguments
argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
} else if (file.exists("./src/utils.r")) {
  source("./src/utils.r", chdir = FALSE)
}

# Default arguments for interactive testing
if (interactive()) {
  source("./src/utils.r", chdir = FALSE)
  argvs$distances <- "data/TypeToTypeDistances.csv"
  argvs$annotations <- "data/visual_neurons_anno.csv"
  argvs$output <- "similarity_tree.pdf"
}

# Validate required arguments
if (is.null(argvs$distances)) stop("--distances argument is required")
if (is.null(argvs$annotations)) stop("--annotations argument is required")
if (is.null(argvs$output)) argvs$output <- "similarity_tree.pdf"

# Load data
cat("Loading distance matrix...\n")
sim_o <- read.csv(argvs$distances, row.names = 1)

cat("Loading annotations...\n")
opc <- read.csv(argvs$annotations)
opc <- subset(opc, temporal_id != 0)

# Define cell types to remove (from original script)
to_remove <- c("C2", "C3", "T2", "T2a", "T3", "T4a", "T4b", "T4c", "T4d",
               "T5a", "T5b", "T5c", "T5d", "L1", "L2", "L3", "L4", "L5",
               "Y1", "Y3", "Y4", "Y11", "Y12", "Lawf1", "Lawf2", "LPi01",
               "LPi02", "LPi03", "LPi04", "LPi05", "LPi06", "LPi07",
               "LPi08", "LPi10", "LPi11", "LPi12", "LPi13", "LPi14", "LPi15",
               paste0("Li0", c(1:5, 7:9)), paste0("Li", c(10:13, 15:33)), "CT1",
               "Am1", "LLPt", "PDt", "LMa1", "LMa5", paste0("LMt", 1:4),
               "Pm12", "Pm13", "Pm14", "Sm38", "Sm39", "Tlp1", "Tlp4", "Tlp5",
               "Tlp14")

# Filter distance matrix
to_keep_r <- setdiff(row.names(sim_o), to_remove)
to_keep_c <- setdiff(colnames(sim_o), to_remove)
sim <- sim_o[to_keep_r, to_keep_c] %>% dist()

cat("Generating hierarchical clustering...\n")

# Create lookup tables
lut <- opc$temporal_label
names(lut) <- opc$cell_type

dlut <- opc$temporal_id
names(dlut) <- opc$cell_type

func_lut <- opc$func
names(func_lut) <- opc$cell_type

# Extract color palettes from existing scale functions
temporal_scale <- scale_color_nk2023()
temporal_colors <- temporal_scale$palette.cache
if (is.null(temporal_colors)) {
  # Extract from scale_color_manual values
  temporal_colors <- c("Hth" = "#C76F6B",
                      "Between Hth & Hth/Opa" = "#EFA900",
                      "Hth/Opa" = "#EFC8B9",
                      "Opa/Erm" = "#FEE699",
                      "Erm/Ey" = "#79A68C",
                      "Ey/Hbn" = "#82C9C5",
                      "Hbn/Opa/Slp" = "#B4C6E6",
                      "Slp/D" = "#A6A1CD",
                      "D/BH-1" = "#A46690")
}

functional_scale <- scale_color_subsystem()
functional_colors <- functional_scale$palette.cache
if (is.null(functional_colors)) {
  # Extract from scale_color_manual values
  functional_colors <- c(
    Color = "#7FC97F",
    Object = "#BEAED4",
    Motion = "#FDC086",
    Luminance = "#FFFF99",
    Unannotated = "#386CB0",
    Polarization = "#F0027F",
    Form = "#BF5B17"
  )
}

cat("Creating temporal origin tree...\n")
# Tree 1: Temporal origin
ttree <- hclust(sim) %>%
  as.phylo() %>%
  as_tibble() %>%
  mutate(
    tlabel = ifelse(label %in% names(lut), lut[label], "Others")
  ) %>%
  mutate(
    tlabel = factor(tlabel, levels = c("Others", names(temporal_colors)))
  ) %>%
  as.treedata() %>%
  ggtree(layout='circular', color = "black") +
  geom_point(aes(color = tlabel), size = 0) +
  geom_tiplab(aes(color = tlabel), show.legend = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = c(temporal_colors, "Others" = "black"), 
                     breaks = names(temporal_colors)) +
  theme(
    plot.background = element_rect(fill = "#ffffff", color = "#ffffff"),
    panel.background = element_rect(fill = "#ffffff", color = "#ffffff"),
    legend.background = element_rect(fill = "#ffffff"),
    legend.text = element_text(color = "black"),
    legend.title = element_blank()
  )

cat("Creating functional subsystem tree...\n")
# Tree 2: Functional subsystem
func_tree <- hclust(sim) %>%
  as.phylo() %>%
  as_tibble() %>%
  mutate(
    functional = ifelse(label %in% names(func_lut), func_lut[label], "Others")
  ) %>%
  mutate(
    functional = factor(functional, levels = c("Others", names(functional_colors)))
  ) %>%
  as.treedata() %>%
  ggtree(layout='circular', color = "black") +
  geom_point(aes(color = functional), size = 0) +
  geom_tiplab(aes(color = functional), show.legend = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = c(functional_colors, "Others" = "black"), na.translate = FALSE) +
  theme(
    plot.background = element_rect(fill = "#ffffff", color = "#ffffff"),
    panel.background = element_rect(fill = "#ffffff", color = "#ffffff"),
    legend.background = element_rect(fill = "#ffffff"),
    legend.text = element_text(color = "black"),
    legend.title = element_blank()
  )

cat("Combining plots and saving...\n")
# Combine both plots
combined_plot <- ttree / func_tree

# Save the plot
ggsave(
  argvs$output,
  combined_plot,
  device = cairo_pdf,
  width = 12, 
  height = 14
)

cat(sprintf("Similarity tree saved to: %s\n", argvs$output))