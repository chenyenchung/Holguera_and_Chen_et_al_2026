#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

# Interactive mode for development
if (interactive()) {
  source("./src/utils.r", chdir = FALSE)
  argvs$connections <- 'data/connections_princeton_no_threshold.csv.gz'
  argvs$celltypes <- 'data/consolidated_cell_types.csv.gz'
  argvs$rotated_dir <- 'int/idv_mat'
  argvs$viz_meta <- 'data/viz_meta.csv'
  argvs$lc_lplc_map <- 'data/lc_lplc_dev_origin.csv'
  argvs$temporal_map <- 'data/temporal_to_dev_origin.csv'
  argvs$visual_anno <- 'data/visual_neurons_anno.csv'
  argvs$output_dir <- 'plots/lc_lplc_dev_origin/version_A'
  argvs$version <- 'A'  # or 'B'
  argvs$syn_threshold <- 5
} else {
  # Ensure utils.r is sourced
  if (!file.exists("./utils.r") && file.exists("src/utils.r")) {
    source("src/utils.r", chdir = FALSE)
  }
}

# Set default parameters
if (is.null(argvs$syn_threshold)) {
  argvs$syn_threshold <- 5
} else {
  argvs$syn_threshold <- as.numeric(argvs$syn_threshold)
}

if (is.null(argvs$version)) {
  argvs$version <- 'A'
}

# Validate input files
if (is.null(argvs$connections) || !file.exists(argvs$connections)) {
  stop("Connection file not found: ", argvs$connections)
}

if (is.null(argvs$celltypes) || !file.exists(argvs$celltypes)) {
  stop("Cell type annotation file not found: ", argvs$celltypes)
}

if (is.null(argvs$viz_meta) || !file.exists(argvs$viz_meta)) {
  stop("Visualization metadata file not found: ", argvs$viz_meta)
}

if (is.null(argvs$rotated_dir) || !dir.exists(argvs$rotated_dir)) {
  stop("Rotated matrix directory not found: ", argvs$rotated_dir)
}

if (is.null(argvs$lc_lplc_map) || !file.exists(argvs$lc_lplc_map)) {
  stop("LC/LPLC developmental origin map not found: ", argvs$lc_lplc_map)
}

if (argvs$version == "B") {
  if (is.null(argvs$temporal_map) || !file.exists(argvs$temporal_map)) {
    stop("Temporal mapping file required for Version B: ", argvs$temporal_map)
  }
  if (is.null(argvs$visual_anno) || !file.exists(argvs$visual_anno)) {
    stop("Visual neuron annotation file required for Version B: ", argvs$visual_anno)
  }
}

if (is.null(argvs$output_dir)) {
  argvs$output_dir <- sprintf('plots/lc_lplc_dev_origin/version_%s', argvs$version)
}

# Create output directory if it doesn't exist
if (!dir.exists(argvs$output_dir)) {
  dir.create(argvs$output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", argvs$output_dir))
}

cat("========================================\n")
cat("LC/LPLC Partner Inputs Visualization\n")
cat("Visualizing inputs TO LC/LPLC partners\n")
cat(sprintf("Version: %s\n", argvs$version))
cat("========================================\n\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Load developmental origin data
#' @param lc_lplc_map_file Path to LC/LPLC developmental origin mapping
#' @param temporal_map_file Path to temporal-to-dev-origin mapping (optional)
#' @param visual_anno_file Path to visual neuron annotations (optional)
#' @return List with lc_lplc_map, temporal_map, and visual_anno
load_dev_origin_data <- function(lc_lplc_map_file, temporal_map_file = NULL, visual_anno_file = NULL) {
  lc_lplc_map <- fread(lc_lplc_map_file)
  result <- list(lc_lplc_map = lc_lplc_map)

  if (!is.null(temporal_map_file)) {
    result$temporal_map <- fread(temporal_map_file)
  }
  if (!is.null(visual_anno_file)) {
    result$visual_anno <- fread(visual_anno_file)
  }

  return(result)
}

#' Identify partner neurons (root set A) and annotate with color variable
#' Root set A: neurons that synapse onto LC/LPLCs
#' @param connections Connection data.table
#' @param target_lc_lplcs Vector of target LC/LPLC types
#' @param syn_threshold Minimum synapse count
#' @param coord Coordinate data.table
#' @param lc_lplc_map LC/LPLC developmental origin mapping
#' @param visual_anno Visual neuron annotations (for Version B)
#' @param temporal_map Temporal-to-dev-origin mapping (for Version B)
#' @param version "A" or "B"
#' @param color_by "origin" or "lc_lplc_type"
#' @return Data.table with post_root_id (from synapse set B) and color_var columns
identify_and_annotate_partners <- function(connections, target_lc_lplcs, syn_threshold, coord,
                                            lc_lplc_map, visual_anno = NULL, temporal_map = NULL,
                                            version = "A", color_by = "origin") {
  # Step 1: Identify root set A (partners of LC/LPLCs)
  # These are neurons that synapse ONTO target LC/LPLCs
  ol_npils <- c("ME_L", "ME_R", "LO_L", "LO_R", "LOP_L", "LOP_R")
  partner_dt <- connections[
    post_type %in% target_lc_lplcs & syn_count > syn_threshold & neuropil %in% ol_npils,
    .(pre_root_id, post_type)
  ]

  if (nrow(partner_dt) == 0) {
    return(data.table(post_root_id = character(0), color_var = character(0)))
  }

  # Merge LC/LPLC developmental origins
  lc_lplc_origins <- merge(
    data.table(lc_lplc_type = target_lc_lplcs),
    lc_lplc_map,
    by.x = "lc_lplc_type",
    by.y = "lc_lplc_type",
    all.x = TRUE
  )
  lc_lplc_origins[is.na(dev_origin), dev_origin := "unknown"]

  if (version == "A") {
    if (color_by == "origin") {
      # Color by LC/LPLC developmental origin that the partner connects to
      partner_dt <- merge(partner_dt, lc_lplc_origins[, .(lc_lplc_type, dev_origin)],
                         by.x = "post_type", by.y = "lc_lplc_type")
      # For multi-origin partners, use majority
      partner_colors <- partner_dt[, .(color_var = names(sort(table(dev_origin), decreasing=TRUE))[1]),
                                   by = pre_root_id]
    } else if (color_by == "lc_lplc_type") {
      # Color by specific LC/LPLC type that the partner connects to
      partner_colors <- partner_dt[, .(color_var = names(sort(table(post_type), decreasing=TRUE))[1]),
                                   by = pre_root_id]
    }
    # Rename for synapse set B (where partners are postsynaptic)
    setnames(partner_colors, "pre_root_id", "post_root_id")

  } else if (version == "B") {
    # Color by temporal origin of the partner neuron itself
    # Get partner types from coordinate data
    partner_types <- unique(coord[post_root_id %in% unique(partner_dt$pre_root_id), .(post_root_id, post_type)])
    partner_annotated <- merge(partner_types, visual_anno[, .(cell_type, temporal_label)],
                              by.x = "post_type", by.y = "cell_type", all.x = TRUE)
    partner_annotated[is.na(temporal_label), temporal_label := "unknown"]

    partner_annotated <- merge(partner_annotated, temporal_map, by = "temporal_label", all.x = TRUE)
    partner_annotated[is.na(dev_origin_upstream), dev_origin_upstream := "unknown_temporal"]

    partner_colors <- partner_annotated[, .(post_root_id, color_var = dev_origin_upstream)]
  }

  return(partner_colors)
}

#' Annotate and stratified subsample coordinates
#' @param coord Coordinate data.table (synapse set B)
#' @param partner_colors Partner neuron color assignments (by post_root_id)
#' @param sample_size Total synapses to sample
#' @return Subsampled data.table with color_var column
annotate_and_subsample <- function(coord, partner_colors, sample_size) {
  # Merge by post_root_id since we're looking at inputs TO partners
  coord <- merge(coord, partner_colors, by = "post_root_id", all.x = TRUE)
  coord[is.na(color_var), color_var := "non_partner"]

  # Stratified sampling
  color_counts <- coord[, .N, by = color_var]
  n_total <- nrow(coord)
  if (n_total <= sample_size) return(coord)

  color_counts[, target_n := round(sample_size * N / n_total)]
  # Adjust total
  diff <- sample_size - sum(color_counts$target_n)
  if (diff != 0) color_counts[which.max(N), target_n := target_n + diff]

  set.seed(1)
  sampled <- rbindlist(lapply(color_counts$color_var, function(cv) {
    dt <- coord[color_var == cv]
    n <- color_counts[color_var == cv, target_n]
    if (nrow(dt) <= n) dt else dt[sample(.N, n)]
  }))

  return(sampled)
}

#' Calculate depth density for color groups
#' @param coord_raw Unsubsampled coordinate data with color_var (synapse set B)
#' @param meta Visualization metadata
#' @return Data.table with x (depth), y (density), color_var
calculate_depth_density <- function(coord_raw, meta) {
  # Use postsynaptic depth (where partners receive inputs)
  depth_data <- coord_raw[, .(color_var, depth = post_rz)]
  if (meta$zinvert) depth_data[, depth := depth * -1]

  depth_by_color <- split(depth_data, depth_data$color_var, drop = TRUE)
  den_grid <- seq(meta$minz, meta$maxz, length.out = 1024)

  interpolated <- rbindlist(lapply(names(depth_by_color), function(cv) {
    d <- density(depth_by_color[[cv]]$depth)
    y_interp <- approx(d$x, d$y, xout = den_grid, rule = 2)$y
    data.table(x = den_grid, y = y_interp, color_var = cv)
  }))

  return(interpolated)
}

#' Create combined scatter + density plot
#' @param coord_sub Subsampled coordinates with color_var (synapse set B)
#' @param coord_raw Unsubsampled coordinates for density (synapse set B)
#' @param meta Visualization metadata
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param color_scale ggplot color scale function
#' @return List with plot and legend
create_combined_plot <- function(coord_sub, coord_raw, meta, title, subtitle, color_scale) {
  # Use postsynaptic coordinates (where partners receive inputs)
  x_axis <- paste0("post_", meta$x_axis)
  y_axis <- paste0("post_", meta$y_axis)

  # Scatter plot - ensure non_partner points are plotted first (underneath)
  coord_sub[, plot_order := ifelse(color_var == "non_partner", 1, 2)]
  setorder(coord_sub, plot_order, color_var)
  scatter <- ggplot(coord_sub, aes(x = .data[[x_axis]], y = .data[[y_axis]], color = color_var)) +
    rasterise(geom_point(size = 0.5, alpha = 0.6), dpi = 450) +
    color_scale +
    labs(title = title, subtitle = subtitle) +
    theme_ih2025()

  # Apply axis limits with reversals
  scale_axis_1 <- get(meta$axis_1_func)
  scale_axis_2 <- get(meta$axis_2_func)
  min1 <- if (grepl("reverse", meta$axis_1_func)) meta$max1 else meta$min1
  max1 <- if (grepl("reverse", meta$axis_1_func)) meta$min1 else meta$max1
  min2 <- if (grepl("reverse", meta$axis_2_func)) meta$max2 else meta$min2
  max2 <- if (grepl("reverse", meta$axis_2_func)) meta$min2 else meta$max2

  scatter <- scatter + scale_axis_1(limits = c(min1, max1)) + scale_axis_2(limits = c(min2, max2))

  # Extract and remove legend
  legend <- get_plot_component(scatter, "guide-box", return_all = TRUE)
  if (!inherits(legend, "grob")) {
    legend_list <- legend[sapply(legend, function(x) !"zeroGrob" %in% class(x))]
    if (length(legend_list) > 0) {
      legend <- legend_list[[1]]
    }
  }
  scatter <- scatter + theme(legend.position = "none")

  # Density plot
  density_data <- calculate_depth_density(coord_raw[color_var != "non_partner"], meta)
  density_plot <- ggplot(density_data, aes(x = x, y = y, color = color_var)) +
    geom_line(linewidth = 1) +
    theme_ih2025() +
    color_scale +
    guides(color = "none")

  if (meta$zid == "y") density_plot <- density_plot + coord_flip()

  # Combine
  if (meta$outlayout == "landscape") {
    combined <- (scatter | density_plot) + plot_layout(widths = c(meta$outr1, meta$outr2))
  } else {
    combined <- (scatter / density_plot) + plot_layout(heights = c(meta$outr1, meta$outr2))
  }

  return(list(plot = combined, legend = legend))
}

# ============================================================================
# MAIN SCRIPT
# ============================================================================

# Load data
cat("Loading data...\n")
connections <- fread(argvs$connections, colClasses = c(pre_root_id = "character", post_root_id = "character"))
cat(sprintf("  Loaded connections: %d rows\n", nrow(connections)))

type_ann <- fread(argvs$celltypes, colClasses = c(root_id = "character"))
cat(sprintf("  Loaded cell type annotations: %d neurons\n", nrow(type_ann)))

viz_meta <- fread(argvs$viz_meta)
cat(sprintf("  Loaded visualization metadata: %d neuropils\n", nrow(viz_meta)))

# Load developmental origin mappings
cat(sprintf("Loading developmental origin data (Version %s)...\n", argvs$version))
dev_maps <- load_dev_origin_data(
  argvs$lc_lplc_map,
  if (argvs$version == "B") argvs$temporal_map else NULL,
  if (argvs$version == "B") argvs$visual_anno else NULL
)
cat(sprintf("  LC/LPLC origin mappings: %d types\n", nrow(dev_maps$lc_lplc_map)))

# Normalize root IDs by stripping prefix
cat("\nNormalizing root IDs...\n")
connections[, pre_root_id := sub("^720575940", "", pre_root_id)]
connections[, post_root_id := sub("^720575940", "", post_root_id)]
type_ann[, root_id := sub("^720575940", "", root_id)]

# Annotate connections with cell types
cat("Annotating connections with cell types...\n")

# Merge presynaptic types
connections <- merge(connections, type_ann[, .(root_id, primary_type)],
                     by.x = "pre_root_id", by.y = "root_id", all.x = TRUE)
setnames(connections, "primary_type", "pre_type")

# Merge postsynaptic types
connections <- merge(connections, type_ann[, .(root_id, primary_type)],
                     by.x = "post_root_id", by.y = "root_id", all.x = TRUE)
setnames(connections, "primary_type", "post_type")

# Filter for LC/LPLC postsynaptic types
connections <- connections[grepl("^LC[0-9]+|^LPLC[0-9]+", post_type)]

cat(sprintf("  Filtered to %d LC/LPLC connections\n", nrow(connections)))

# Get unique LC/LPLC types
unique_lc_lplcs <- sort(unique(connections$post_type))
cat(sprintf("\nFound %d LC/LPLC types total\n", length(unique_lc_lplcs)))

# Get neuropils from metadata
neuropils <- unique(viz_meta$neuropil)
cat(sprintf("Processing %d neuropils: %s\n\n",
            length(neuropils), paste(neuropils, collapse=", ")))

# Define plot groups
cat("Defining plot groups...\n")

# Get LC/LPLC types for each origin
slp_dll_types <- c("LC6", "LC9", "LC10a", "LC11", "LC21")
topc_types <- c("LC13", "LC15", "LC16", "LC24", "LC26")
vtipc_types <- c("LC4", "LC12", "LC17")
cb_types <- c("LPLC1", "LPLC2", "LPLC4", "LPLC22")
unknown_types <- setdiff(unique_lc_lplcs, c(slp_dll_types, topc_types, vtipc_types, cb_types))

plot_groups <- list(
  slp_dll = list(
    name = "Slp_vs_Dll",
    types = slp_dll_types,
    color_by = "origin",
    origins = c("Slp", "Dll")
  ),
  tOPC = list(
    name = "tOPC",
    types = topc_types,
    color_by = "lc_lplc_type",
    origins = "tOPC"
  ),
  vtIPC = list(
    name = "vtIPC",
    types = vtipc_types,
    color_by = "lc_lplc_type",
    origins = "vtIPC"
  ),
  central_brain = list(
    name = "central_brain",
    types = cb_types,
    color_by = "lc_lplc_type",
    origins = "central_brain"
  ),
  unknown = list(
    name = "unknown",
    types = unknown_types,
    color_by = "lc_lplc_type",
    origins = "unknown"
  ),
  mega = list(
    name = "ALL_origins",
    types = unique_lc_lplcs,
    color_by = "origin",
    origins = c("Slp", "Dll", "tOPC", "vtIPC", "central_brain", "unknown")
  )
)

for (gname in names(plot_groups)) {
  cat(sprintf("  %s: %d types\n", plot_groups[[gname]]$name, length(plot_groups[[gname]]$types)))
}
cat("\n")

# Main visualization loop
plots_created <- 0
plots_skipped <- 0

cat("Generating plots...\n")
cat(sprintf("Synapse threshold: >%d synapses\n", argvs$syn_threshold))
cat(sprintf("Version: %s\n\n", argvs$version))

for (group_id in names(plot_groups)) {
  group <- plot_groups[[group_id]]

  cat(sprintf("Processing group: %s (%d types)\n", group$name, length(group$types)))

  # Skip if no types in this group
  if (length(group$types) == 0) {
    cat("  No types in group, skipping\n\n")
    next
  }

  for (np in neuropils) {
    # Load coordinates
    coord_file <- file.path(argvs$rotated_dir, paste0(np, "_rotated.csv.gz"))
    if (!file.exists(coord_file)) {
      cat(sprintf("  %s: coordinate file not found (skipped)\n", np))
      plots_skipped <- plots_skipped + 1
      next
    }

    coord <- fread(coord_file, colClasses = c(pre_root_id = "character", post_root_id = "character"))

    # Normalize root IDs
    coord[, pre_root_id := sub("^720575940", "", pre_root_id)]
    coord[, post_root_id := sub("^720575940", "", post_root_id)]

    # Step 1: Identify root set A (partners of LC/LPLCs) and get their color assignments
    # This returns a mapping of post_root_id -> color_var
    partner_colors <- identify_and_annotate_partners(
      connections, group$types, argvs$syn_threshold, coord,
      dev_maps$lc_lplc_map, dev_maps$visual_anno, dev_maps$temporal_map,
      version = argvs$version, color_by = group$color_by
    )

    # If no partners found, skip
    if (nrow(partner_colors) == 0) {
      cat(sprintf("  %s: no partners for this group (skipped)\n", np))
      plots_skipped <- plots_skipped + 1
      next
    }

    # Step 2: Keep ALL synapses for visualization context
    # Stratified sampling will maintain ratio of synapse set B vs rest
    # Keep copy of all coordinates (unsampled) for density calculation
    coord_raw <- copy(coord)

    # Subsample with stratified sampling (maintains ratio of synapse set B vs rest)
    sample_size <- ifelse(grepl("^ME_", np), 40000, 20000)
    coord_sub <- annotate_and_subsample(coord, partner_colors, sample_size)
    coord_raw <- merge(coord_raw, partner_colors, by = "post_root_id", all.x = TRUE)
    coord_raw[is.na(color_var), color_var := "non_partner"]

    # Check if we have any partner synapses after annotation
    n_partner <- sum(coord_sub$color_var != "non_partner")
    if (n_partner == 0) {
      cat(sprintf("  %s: no inputs to partners (skipped)\n", np))
      plots_skipped <- plots_skipped + 1
      next
    }

    # Get metadata
    meta <- viz_meta[neuropil == np]
    if (nrow(meta) == 0) {
      cat(sprintf("  %s: no metadata (skipped)\n", np))
      plots_skipped <- plots_skipped + 1
      next
    }
    meta <- meta[1, ]

    # Determine color scale
    if (argvs$version == "A") {
      if (group$color_by == "origin") {
        color_scale <- scale_color_dev_origin()
      } else {  # lc_lplc_type
        base_colors <- c(tOPC = "#4DAF4A", vtIPC = "#984EA3",
                        central_brain = "#FF7F00", unknown = "#A65628")
        base_color <- base_colors[group$origins[1]]
        color_scale <- scale_color_lc_lplc_shaded(group$types, base_color)
      }
    } else {  # Version B
      color_scale <- scale_color_temporal_origin()
    }

    # Create plot
    title <- sprintf("Inputs to %s partners in %s", group$name, np)
    subtitle <- sprintf("n = %d synapses | %d to partners (%.1f%%) | Version %s",
                       nrow(coord_sub), n_partner, 100*n_partner/nrow(coord_sub), argvs$version)

    plot_result <- create_combined_plot(coord_sub, coord_raw, meta, title, subtitle, color_scale)

    # Save
    output_file <- file.path(argvs$output_dir,
                             sprintf("%s_%s_version%s.pdf", np, group$name, argvs$version))
    ggsave(output_file, plot_result$plot, width = meta$outd1, height = meta$outd2)

    legend_file <- file.path(argvs$output_dir,
                            sprintf("%s_%s_version%s_legend.pdf", np, group$name, argvs$version))
    ggsave(legend_file, plot_result$legend)

    cat(sprintf("  %s: saved (%d synapses, %.1f%% partners)\n",
                np, nrow(coord_sub), 100*n_partner/nrow(coord_sub)))

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
cat(sprintf("\nExpected total: %d (%d groups × %d neuropils)\n",
            length(plot_groups) * length(neuropils),
            length(plot_groups), length(neuropils)))
