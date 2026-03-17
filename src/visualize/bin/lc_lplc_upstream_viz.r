#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

# Interactive mode for development
if (interactive()) {
  source("./src/utils.r", chdir = FALSE)
  argvs$partner_summary <- 'lc_lplc_presynaptic_partners.csv'
  argvs$anno <- 'data/visual_neurons_anno.csv'
  argvs$rotated_dir <- 'int/idv_mat'
  argvs$output_dir <- 'plots/lc_lplc_upstream'
  argvs$viz_meta <- 'data/viz_meta.csv'
} else {
  # Ensure utils.r is sourced
  if (!file.exists("./utils.r") && file.exists("src/utils.r")) {
    source("src/utils.r", chdir = FALSE)
  }
}

# Validate input files
if (is.null(argvs$partner_summary) || !file.exists(argvs$partner_summary)) {
  stop("Partner summary file not found: ", argvs$partner_summary)
}

if (is.null(argvs$anno) || !file.exists(argvs$anno)) {
  stop("Annotation file not found: ", argvs$anno)
}

if (is.null(argvs$viz_meta) || !file.exists(argvs$viz_meta)) {
  stop("Visualization metadata file not found: ", argvs$viz_meta)
}

if (is.null(argvs$rotated_dir) || !dir.exists(argvs$rotated_dir)) {
  stop("Rotated matrix directory not found: ", argvs$rotated_dir)
}

if (is.null(argvs$output_dir)) {
  argvs$output_dir <- 'plots/lc_lplc_upstream'
}

# Create output directory if it doesn't exist
if (!dir.exists(argvs$output_dir)) {
  dir.create(argvs$output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", argvs$output_dir))
}

cat("========================================\n")
cat("LC/LPLC Upstream Visualization\n")
cat("========================================\n\n")

# Load data
cat("Loading data...\n")
partner_summary <- fread(argvs$partner_summary)
cat(sprintf("  Loaded partner summary: %d rows\n", nrow(partner_summary)))

anno <- fread(argvs$anno)
cat(sprintf("  Loaded annotations: %d cell types\n", nrow(anno)))

viz_meta <- fread(argvs$viz_meta)
cat(sprintf("  Loaded visualization metadata: %d neuropils\n", nrow(viz_meta)))

# Get unique LC/LPLC types
lc_lplc_types <- sort(unique(partner_summary$post_type))
cat(sprintf("\nFound %d LC/LPLC types to visualize:\n", length(lc_lplc_types)))
cat(sprintf("  %s\n", paste(lc_lplc_types, collapse=", ")))

# Get neuropils from metadata
neuropils <- unique(viz_meta$neuropil)
cat(sprintf("\nProcessing %d neuropils: %s\n\n",
            length(neuropils), paste(neuropils, collapse=", ")))

# Function to create contribution plot
create_contribution_plot <- function(coord, meta, lc_lplc_type, neuropil) {
  # Extract axes from metadata
  x_axis <- paste0("pre_", meta$x_axis)
  y_axis <- paste0("pre_", meta$y_axis)

  # Validate axes exist
  if (!x_axis %in% colnames(coord) || !y_axis %in% colnames(coord)) {
    warning(sprintf("Missing axes for %s: %s or %s", neuropil, x_axis, y_axis))
    return(NULL)
  }

  # Create scatter plot
  p <- ggplot(coord, aes_string(x = x_axis, y = y_axis,
                                 color = "contribution_score")) +
    rasterise(
      geom_point(size = 0.5, alpha = 0.6),
      dpi = 450
    ) +
    scale_color_gradient(
      low = "grey92",
      high = "darkblue",
      name = "Contribution Score\n(Connected / Total)",
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0%", "25%", "50%", "75%", "100%")
    ) +
    labs(
      title = sprintf("Presynaptic Partners of %s in %s",
                      lc_lplc_type, neuropil),
      subtitle = sprintf("n = %d synapses from %d cell types",
                        nrow(coord),
                        length(unique(coord$pre_type)))
    ) +
    coord_fixed(ratio = 1) +
    theme_ih2025()

  # Apply axis scaling functions from metadata (handles reversals)
  # Note: axis_1/2 refer to metadata ordering, not plot x/y axes
  scale_axis_1 <- get(meta$axis_1_func)
  scale_axis_2 <- get(meta$axis_2_func)
  if (grepl("reverse", meta$axis_1_func)) {
    min1 <- meta$max1
    max1 <- meta$min1
    meta$min1 <- min1
    meta$max1 <- max1
  }
  if (grepl("reverse", meta$axis_2_func)) {
    min2 <- meta$max2
    max2 <- meta$min2
    meta$min2 <- min2
    meta$max2 <- max2
  }

  p <- p +
    scale_axis_1(limits = c(meta$min1, meta$max1)) +
    scale_axis_2(limits = c(meta$min2, meta$max2))

  return(p)
}

# Main visualization loop
plots_created <- 0
plots_skipped <- 0

cat("Generating plots...\n")
cat("Progress: [LC/LPLC type] × [neuropil]\n\n")

for (lc_lplc_type in lc_lplc_types) {
  cat(sprintf("Processing %s...\n", lc_lplc_type))

  for (np in neuropils) {
    # Load rotated coordinate matrix
    coord_file <- file.path(argvs$rotated_dir,
                            paste0(np, "_rotated.csv.gz"))

    if (!file.exists(coord_file)) {
      warning(sprintf("  Coordinate file not found: %s", coord_file))
      plots_skipped <- plots_skipped + 1
      next
    }

    coord <- fread(coord_file)

    # Apply filter
    filtered_coord <- filter_lc_lplc_upstream(
      coord, anno, partner_summary, lc_lplc_type, syn_type = "pre"
    )

    # Skip if no data
    if (nrow(filtered_coord) == 0) {
      cat(sprintf("  %s: no data (skipped)\n", np))
      plots_skipped <- plots_skipped + 1
      next
    }

    # Subsample for rendering efficiency
    # Use larger sample for ME (Medulla), smaller for others
    sample_size <- if (grepl("^ME_", np)) 40000 else 20000
    n_before <- nrow(filtered_coord)
    filtered_coord <- rsubsample(filtered_coord, n = sample_size, seed = 1, syn_type = "pre")
    n_after <- nrow(filtered_coord)

    # Get neuropil metadata
    meta <- viz_meta[neuropil == np]

     if (nrow(meta) == 0) {
      warning(sprintf("  No metadata for %s (skipped)", np))
      plots_skipped <- plots_skipped + 1
      next
    }

    # Extract first row (should only be one)
    meta <- meta[1, ]

    # Create plot
    plot <- create_contribution_plot(
      filtered_coord, meta, lc_lplc_type, np
    )

    if (is.null(plot)) {
      plots_skipped <- plots_skipped + 1
      next
    }

    # Save plot
    output_file <- file.path(
      argvs$output_dir,
      sprintf("%s_%s_upstream_contribution.pdf", np, lc_lplc_type)
    )

    # Use dimensions from metadata (landscape: 10x6, portrait: 5x10)
    ggsave(output_file, plot, width = meta$outd1, height = meta$outd2)
    if (n_before > n_after) {
      cat(sprintf("  %s: saved (%d→%d synapses, subsampled)\n", np, n_before, n_after))
    } else {
      cat(sprintf("  %s: saved (%d synapses)\n", np, n_after))
    }
    plots_created <- plots_created + 1
  }

  cat("\n")
}

# Summary
cat("========================================\n")
cat("Visualization complete!\n")
cat("========================================\n")
cat(sprintf("Plots created: %d\n", plots_created))
cat(sprintf("Plots skipped (no data): %d\n", plots_skipped))
cat(sprintf("Output directory: %s\n", argvs$output_dir))
cat(sprintf("\nTotal files expected: %d (48 LC/LPLC × 6 neuropils)\n",
            length(lc_lplc_types) * length(neuropils)))
