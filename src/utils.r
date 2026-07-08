# Consolidated utility functions for flyem analysis
# This file contains shared utility functions used across multiple R scripts

# ============================================================================
# PLOTTING UTILITIES
# ============================================================================

#' Custom ggplot theme for flyem visualizations
theme_ih2025 <- function(xlim = NULL, ylim = NULL) {
  f <- theme_minimal() %+replace%
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(face = "italic"),
      axis.text = element_blank(),
      axis.title = element_blank(),
      strip.background = element_rect(color = "black", fill = "transparent"),
      strip.text = element_text(size = 14, face = "bold"),
      panel.background = element_rect(color = "black", fill = "transparent"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      legend.byrow = TRUE
    )
  return(f)
}

# ============================================================================
# COLOR SCALE FUNCTIONS
# ============================================================================

#' Color scale for NK 2023 data
scale_color_nk2023 <- function() {
  f <- scale_color_manual(
    values = c("Hth" = "#C76F6B",
               "Between Hth & Hth/Opa" = "#EFA900",
               "Hth/Opa" = "#EFC8B9",
               "Opa/Erm" = "#FEE699",
               "Erm/Ey" = "#79A68C",
               "Ey/Hbn" = "#82C9C5",
               "Hbn/Opa/Slp" = "#B4C6E6",
               "Slp/D" = "#A6A1CD",
               "D/BH-1" = "#A46690")
  )
  return(f)
}

ih2025_colors <- function() {
  c(
    "Early" = "#8DC63F",
    "Late" = "#7F55A2"
  )
}

scale_color_ih2025 <- function() {
  f <- scale_color_manual(
    values = ih2025_colors()
  )
  return(f)
}

#' Color scale for subsystem annotations
scale_color_subsystem <- function() {
  f <- scale_color_manual(
    values = c(
      Color = "#7FC97F",
      Object = "#BEAED4",
      Motion = "#FDC086",
      Luminance = "#FFFF99",
      Unannotated = "#386CB0",
      Polarization = "#F0027F",
      Form = "#BF5B17"
    )
  )
  return(f)
}

#' Fill scale for subsystem annotations
scale_fill_subsystem <- function() {
  f <- scale_fill_manual(
    values = c(
      Color = "#7FC97F",
      Object = "#BEAED4",
      Motion = "#FDC086",
      Luminance = "#FFFF99",
      Unannotated = "#386CB0",
      Polarization = "#F0027F",
      Form = "#BF5B17"
    )
  )
  return(f)
}

#' Cell type colors
type_colors <- function() {
  hsort <- function(x) {
    sort_key <- c(
      "L", "A", "B", "C", "D", "M", "E", "F", "G", "H", "I", "J", "K",
      "N", "O", "P", "Q", "R", "T", "S", "U", "V", "W", "X", "Y", "Z"
    )
    xl <- split(x, substr(x, 1, 1))
    xl <- lapply(xl, sort)
    use_key <- sort_key[sort_key %in% names(xl)]
    out <- unlist(xl[use_key])
    return(out)
  }
  types <- unique(opc_anno$cell_type)
  types <- hsort(types)
  
  hue_pal <- function(
    h = c(0, 360) + 15, c = 100, l = 65, h.start = 0, direction = 1
  ) {
    function(n) {
      if (n == 0) return(character(0))
      hues <- seq(h[1], h[2], length.out = n + 1)[1:n]
      hues <- (hues + h.start) %% 360
      if (direction == -1) hues <- rev(hues)
      grDevices::hcl(h = hues, c = c, l = l)
    }
  }
  colors <- hue_pal()(length(types))
  names(colors) <- types
  colors
}

#' Color scale for cell types
scale_color_type <- function() {
  return(
    scale_color_manual(values = type_colors())
  )
}

#' Color scale for medulla reference groups plus T1
scale_color_medulla_reference_t1 <- function() {
  temporal_colors <- ih2025_colors()
  cell_type_colors <- type_colors()
  scale_color_manual(
    values = c(
      "Other" = "grey92",
      "Dm reference" = temporal_colors[["Late"]],
      "Pm reference" = temporal_colors[["Early"]],
      "T1" = cell_type_colors[["T1"]]
    ),
    breaks = c("T1", "Dm reference", "Pm reference", "Other")
  )
}

#' Color scale for spatial origin annotations
scale_color_spatial_origin <- function() {
  spatial_origin_colors <- c(
    Vsx = "#1B9E77",
    Optix = "#7570B3",
    Dpp = "#D95F02",
    "Vsx/Optix" = "#4B8A95",
    "Vsx/Dpp" = "#7A8B3D",
    "Optix/Dpp" = "#A0685B",
    "Vsx/Optix/Dpp" = "#5F5F5F",
    unknown = "grey80"
  )

  f <- scale_color_manual(
    values = spatial_origin_colors,
    breaks = names(spatial_origin_colors)
  )
  return(f)
}

spatial_origin_color_values <- function() {
  c(
    Vsx = "#1B9E77",
    Optix = "#7570B3",
    Dpp = "#D95F02",
    "Vsx/Optix" = "#4B8A95",
    "Vsx/Dpp" = "#7A8B3D",
    "Optix/Dpp" = "#A0685B",
    "Vsx/Optix/Dpp" = "#5F5F5F",
    unknown = "grey80"
  )
}

scale_color_medulla_hth_spatial <- function() {
  temporal_colors <- ih2025_colors()
  scale_color_manual(
    values = c(
      "Other" = "grey92",
      "Dm reference" = temporal_colors[["Late"]],
      "Pm reference" = temporal_colors[["Early"]],
      spatial_origin_color_values()
    ),
    breaks = c(
      names(spatial_origin_color_values()),
      "Dm reference",
      "Pm reference",
      "Other"
    )
  )
}

#' Color scale for developmental origins (Version A - origin-level)
#' Used for Slp/Dll comparison and mega plot
scale_color_dev_origin <- function() {
  f <- scale_color_manual(
    name = "Developmental Origin",
    values = c(
      "tOPC (Slp)" = "#E41A1C",           # Red
      "tOPC (Dll)" = "#377EB8",           # Blue
      "tOPC" = "#4DAF4A",          # Green
      "vtIPC" = "#984EA3",         # Purple
      "central_brain" = "#FF7F00", # Orange
      "unknown" = "#A65628",       # Brown
      "non_partner" = "grey92"     # Light grey
    )
  )
  return(f)
}

#' Color scale for temporal origins (Version B)
#' Maps upstream neuron temporal windows to developmental categories
scale_color_temporal_origin <- function() {
  f <- scale_color_manual(
    name = "Upstream Temporal Origin",
    values = c(
      "early_temporal" = "#79A68C",   # Green
      "mid_temporal" = "#EFA900",     # Orange
      "late_temporal" = "#7F55A2",    # Purple
      "unknown_temporal" = "#A65628", # Brown
      "non_partner" = "grey92"        # Light grey
    )
  )
  return(f)
}

#' Color scale for LC/LPLC types (Version A - shaded by specific type)
#' Generates gradient shades of a base color for each LC/LPLC type
#' @param lc_lplc_types Vector of LC/LPLC type names
#' @param base_color Base color for the gradient (hex code)
#' @return ggplot scale_color_manual object
scale_color_lc_lplc_shaded <- function(lc_lplc_types, base_color) {
  # Generate shades using colorRampPalette
  n_types <- length(lc_lplc_types)

  # Convert hex to RGB
  base_rgb <- col2rgb(base_color) / 255

  # Create lighter and darker versions
  light_rgb <- pmin(base_rgb * 1.3 + 0.3, 1)
  dark_rgb <- pmax(base_rgb * 0.7, 0)

  light_hex <- rgb(light_rgb[1], light_rgb[2], light_rgb[3])
  dark_hex <- rgb(dark_rgb[1], dark_rgb[2], dark_rgb[3])

  # Generate gradient
  colors <- colorRampPalette(c(light_hex, base_color, dark_hex))(n_types)

  # Create color map
  color_map <- setNames(colors, lc_lplc_types)
  color_map["non_partner"] <- "grey92"

  f <- scale_color_manual(
    name = "LC/LPLC Type",
    values = color_map
  )
  return(f)
}

# ============================================================================
# DATA PROCESSING UTILITIES
# ============================================================================

#' Validate configuration file and check required columns
#' @param file_path Path to the configuration file
#' @param required_cols Vector of required column names
#' @return Data table with validated content
validate_config_file <- function(file_path, required_cols) {
  if (!file.exists(file_path)) stop(sprintf("Config file not found: %s", file_path))
  data <- fread(file_path)
  if (nrow(data) == 0) stop(sprintf("Config file is empty: %s", file_path))
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) stop(sprintf("Missing columns in %s: %s", file_path, paste(missing_cols, collapse=", ")))
  return(data)
}

#' Validate that required functions exist
#' @param preset Preset configuration data
#' @param plot_meta Plot metadata configuration
validate_functions <- function(preset, plot_meta) {
  # Check filter function
  if (!exists(preset$filter_func, mode="function")) 
    stop(sprintf("Filter function not found: %s", preset$filter_func))
  
  # Check color function  
  if (!exists(preset$palette, mode="function"))
    stop(sprintf("Color function not found: %s", preset$palette))
    
  # Check axis functions
  if (!exists(plot_meta$axis_1_func, mode="function"))
    stop(sprintf("Axis function not found: %s", plot_meta$axis_1_func))
  if (!exists(plot_meta$axis_2_func, mode="function"))
    stop(sprintf("Axis function not found: %s", plot_meta$axis_2_func))
}

#' Random subsample data with optional seed
#' @param x Data frame or data.table to subsample
#' @param n Maximum number of rows to keep
#' @param seed Random seed for reproducibility
#' @return Subsampled data
rsubsample <- function(x, n = 1e4, seed = 1, syn_type) {
  type_col <- paste(syn_type, "type", sep = "_")
  all_types <- unique(x[, get(type_col)])
  if (!is.null(nrow(x)) && nrow(x) > n) {
    set.seed(seed)
    if (is.data.table(x)) {
      out <- x[sample(.N, n)]
    } else {
      out <- x[sample(seq_len(nrow(x)), n), , drop = FALSE]
    }
    out_types <- unique(out[, get(type_col)])
    # If rare types were dropped, we re-sample 10 synapses for representation
    # for each type.
    dropped_types <- setdiff(all_types, out_types)
    if (length(dropped_types) > 0) {
      dropped_coords <- x[get(type_col) %in% dropped_types]
      reinsert_coords <- x[, .SD[sample(.N, min(.N, 5))], by = type_col]
      out <- rbind(out, reinsert_coords)
    }
    return(out)
  }
  return(x)
}

# ============================================================================
# FILTERING FUNCTIONS
# ============================================================================

#' Filter for type sparseness and the lack of variation
filter_type <- function(x, syn_type, sparse_limit) {
  x <- x[, .row_id := .I]
  # Filter out types with too few synapses for consistent visualization
  type_col <- paste0(syn_type, "_type")
  type_counts <- x[, .N, by = c(type_col)]
  keep_types <- type_counts[N >= sparse_limit, get(type_col)]
  dropped_types <- type_counts[N < sparse_limit]
  if (nrow(dropped_types) > 0) {
    message(sprintf("Dropped %d %s types with < %d synapses", 
                    nrow(dropped_types), syn_type, sparse_limit))
  }
  x <- x[get(type_col) %in% keep_types]
  
  # Filter out types with insufficient variation for density calculation
  depth_col <- paste0(syn_type, "_rz")
  if (depth_col %in% colnames(x)) {
    depth_variation <- x[, .(unique_depths = length(unique(get(depth_col)))), by = c(type_col)]
    valid_types <- depth_variation[unique_depths >= 2, get(type_col)]
    if (length(valid_types) < nrow(depth_variation)) {
      message(sprintf("Dropped %d types with insufficient depth variation for density calculation", 
                      nrow(depth_variation) - length(valid_types)))
    }
    x <- x[get(type_col) %in% valid_types]
    setorder(x, .row_id)
    x[, .row_id := NULL]
    return(x)
  }
}

#' Filter coordinates by temporal annotations
filter_temporal <- function(coord, ann, syn_type = "pre") {
  by_x <- ifelse(syn_type == "pre", "pre_type", "post_type")
  ann <- ann[temporal_label != "unknown" & Confident_annotation == "Y"]
  ann$temporal_label <- factor(
    ann$temporal_label,
    levels = unique(ann$temporal_label),
    labels = unique(ann$temporal_label)
  )
  coord[, .row_id := .I]
  coord <- merge(
    coord,
    ann[, .(cell_type, temporal_label, Notch, newly_ann, ntype)],
    by.x = by_x,
    by.y = "cell_type"
  )
  setorder(coord, .row_id)
  coord[, .row_id := NULL]
  return(coord)
}



#' Filter coordinates by subsystem annotations
filter_subsystem <- function(coord, ann, syn_type = "pre") {
  by_x <- ifelse(syn_type == "pre", "pre_type", "post_type")
  ann <- ann[func != "unknown" & Confident_annotation == "Y"]
  ann$func <- factor(ann$func)
  coord[, .row_id := .I]
  coord <- merge(
    coord,
    ann[, .(cell_type, func, Notch, newly_ann, ntype)],
    by.x = by_x,
    by.y = "cell_type"
  )
  setorder(coord, .row_id)
  coord[, .row_id := NULL]
  return(coord)
}

#' Filter coordinates by broad temporal annotations
filter_broad <- function(coord, ann, syn_type = "pre") {
  by_x <- ifelse(syn_type == "pre", "pre_type", "post_type")
  ann <- ann[broad_temp %in% c("Early", "Late") & Confident_annotation == "Y"]
  ann$broad_temp <- factor(ann$broad_temp)
  coord[, .row_id := .I]
  coord <- merge(
    coord,
    ann[, .(cell_type, broad_temp, Notch, newly_ann, ntype)],
    by.x = by_x,
    by.y = "cell_type"
  )
  setorder(coord, .row_id)
  coord[, .row_id := NULL]
  return(coord)
}

#' Filter coordinates to one exact cell type
filter_cell_type <- function(coord, ann, syn_type = "pre", cell_type) {
  type_col <- paste0(syn_type, "_type")
  coord <- coord[get(type_col) == cell_type]
  coord[, .row_id := .I]
  coord <- merge(
    coord,
    ann[, .(cell_type)],
    by.x = type_col,
    by.y = "cell_type"
  )
  setorder(coord, .row_id)
  coord[, .row_id := NULL]
  return(coord)
}

annotate_medulla_context <- function(coord, syn_type, target_flag,
                                     target_group) {
  type_col <- paste0(syn_type, "_type")
  is_dm_ref <- grepl("^Dm", coord[[type_col]])
  is_pm_ref <- grepl("^Pm", coord[[type_col]])
  coord[, plot_alpha := fifelse(target_flag, 1, 0.05)]
  coord[, plot_order := fifelse(
    target_flag,
    2L,
    fifelse(is_dm_ref | is_pm_ref, 1L, 0L)
  )]
  coord[, reference_group := fifelse(
    target_flag,
    target_group,
    fifelse(
      is_dm_ref,
      "Dm reference",
      fifelse(is_pm_ref, "Pm reference", "Other")
    )
  )]
  coord
}

#' Annotate coordinates for T1 with medulla reference context
filter_t1 <- function(coord, ann, syn_type = "pre") {
  type_col <- paste0(syn_type, "_type")
  coord <- annotate_medulla_context(
    coord,
    syn_type = syn_type,
    target_flag = coord[[type_col]] == "T1",
    target_group = "T1"
  )
  coord[, reference_group := factor(
    reference_group,
    levels = c("Other", "Dm reference", "Pm reference", "T1")
  )]
  coord[, .row_id := .I]
  coord <- merge(
    coord,
    ann[, .(cell_type)],
    by.x = type_col,
    by.y = "cell_type",
    all.x = TRUE
  )
  setorder(coord, .row_id)
  coord[, .row_id := NULL]
  return(coord)
}

add_spatial_origin <- function(ann) {
  spatial_families <- list(
    Vsx = c("vVsx", "dVsx"),
    Optix = c("vOptix", "dOptix"),
    Dpp = c("vDpp", "dDpp")
  )
  spatial_cols <- unlist(spatial_families, use.names = FALSE)
  spatial_levels <- c(
    "Vsx", "Optix", "Dpp",
    "Vsx/Optix", "Vsx/Dpp", "Optix/Dpp",
    "Vsx/Optix/Dpp", "unknown"
  )
  missing_cols <- setdiff(spatial_cols, colnames(ann))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing spatial origin columns in annotation data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  is_spatial_origin <- function(x) {
    x_chr <- toupper(trimws(as.character(x)))
    !is.na(x_chr) & x_chr %in% c("1", "TRUE", "T", "Y", "YES")
  }

  ann_spatial <- copy(ann)
  for (col in spatial_cols) {
    ann_spatial[, (col) := is_spatial_origin(get(col))]
  }

  family_order <- names(spatial_families)
  family_flag_cols <- paste0(".spatial_", family_order)
  names(family_flag_cols) <- family_order
  for (family in family_order) {
    ann_spatial[
      ,
      (family_flag_cols[[family]]) := rowSums(.SD) > 0,
      .SDcols = spatial_families[[family]]
    ]
  }

  ann_spatial[, spatial_origin := apply(.SD, 1, function(row) {
    origins <- family_order[as.logical(row)]
    if (length(origins) == 0) {
      return("unknown")
    }
    paste(origins, collapse = "/")
  }), .SDcols = family_flag_cols]
  ann_spatial[, spatial_origin := factor(
    spatial_origin,
    levels = spatial_levels
  )]
  ann_spatial
}

#' Annotate exact-Hth cells with medulla reference context
filter_hth_spatial <- function(coord, ann, syn_type = "pre") {
  type_col <- paste0(syn_type, "_type")
  ann_spatial <- add_spatial_origin(ann)
  ann_spatial[, is_hth := !is.na(temporal_label) &
    temporal_label == "Hth" &
    Confident_annotation == "Y"]
  coord[, .row_id := .I]
  coord <- merge(
    coord,
    ann_spatial[, .(cell_type, spatial_origin, is_hth)],
    by.x = type_col,
    by.y = "cell_type",
    all.x = TRUE
  )
  setorder(coord, .row_id)
  coord[, .row_id := NULL]
  coord[is.na(is_hth), is_hth := FALSE]
  coord <- annotate_medulla_context(
    coord,
    syn_type = syn_type,
    target_flag = coord$is_hth,
    target_group = as.character(coord$spatial_origin)
  )
  coord[!is_hth & reference_group %in% c("Dm reference", "Pm reference"),
        spatial_origin := NA]
  coord[, reference_group := factor(
    reference_group,
    levels = c(
      "Other",
      "Dm reference",
      "Pm reference",
      names(spatial_origin_color_values())
    )
  )]
  coord[, is_hth := NULL]
  return(coord)
}


#' Filter coordinates by spatial origin annotations
filter_spatial_origin <- function(coord, ann, syn_type = "pre") {
  by_x <- ifelse(syn_type == "pre", "pre_type", "post_type")
  spatial_families <- list(
    Vsx = c("vVsx", "dVsx"),
    Optix = c("vOptix", "dOptix"),
    Dpp = c("vDpp", "dDpp")
  )
  spatial_cols <- unlist(spatial_families, use.names = FALSE)
  spatial_levels <- c(
    "Vsx", "Optix", "Dpp",
    "Vsx/Optix", "Vsx/Dpp", "Optix/Dpp",
    "Vsx/Optix/Dpp", "unknown"
  )
  missing_cols <- setdiff(spatial_cols, colnames(ann))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing spatial origin columns in annotation data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  is_spatial_origin <- function(x) {
    x_chr <- toupper(trimws(as.character(x)))
    !is.na(x_chr) & x_chr %in% c("1", "TRUE", "T", "Y", "YES")
  }

  ann_spatial <- copy(ann)
  for (col in spatial_cols) {
    ann_spatial[, (col) := is_spatial_origin(get(col))]
  }

  family_order <- names(spatial_families)
  family_flag_cols <- paste0(".spatial_", family_order)
  names(family_flag_cols) <- family_order
  for (family in family_order) {
    ann_spatial[
      ,
      (family_flag_cols[[family]]) := rowSums(.SD) > 0,
      .SDcols = spatial_families[[family]]
    ]
  }

  ann_spatial[, spatial_origin := apply(.SD, 1, function(row) {
    origins <- family_order[as.logical(row)]
    if (length(origins) == 0) {
      return("unknown")
    }
    paste(origins, collapse = "/")
  }), .SDcols = family_flag_cols]
  ann_spatial[, spatial_origin := factor(
    spatial_origin,
    levels = spatial_levels
  )]

  coord[, .row_id := .I]
  coord <- merge(
    coord,
    ann_spatial[, .(cell_type, spatial_origin, Notch, newly_ann, ntype)],
    by.x = by_x,
    by.y = "cell_type"
  )
  setorder(coord, .row_id)
  coord[, .row_id := NULL]
  return(coord)
}

#' Filter coordinates by putative annotations
filter_putative <- function(coord, ann, syn_type = "pre") {
  by_x <- ifelse(syn_type == "pre", "pre_type", "post_type")
  ann <- ann[putative_OPC == TRUE]
  putative_hl_cols <- grep("^putative_hl[0-9]+$", names(ann), value = TRUE)
  merge_cols <- c(
    "cell_type", "Notch", "newly_ann", "ntype", "func", putative_hl_cols
  )
  coord[, .row_id := .I]
  coord <- merge(
    coord,
    ann[, ..merge_cols],
    by.x = by_x,
    by.y = "cell_type"
  )
  setorder(coord, .row_id)
  coord[, .row_id := NULL]
  return(coord)
}

#' Filter coordinates by temporal annotations (all annotated)
filter_temporal_all <- function(coord, ann, syn_type = "pre") {
  by_x <- ifelse(syn_type == "pre", "pre_type", "post_type")
  ann <- ann[Confident_annotation == "Y"]
  ann$temporal_label <- factor(
    ann$temporal_label,
    levels = unique(ann$temporal_label),
    labels = unique(ann$temporal_label)
  )
  coord[, .row_id := .I]
  coord <- merge(
    coord,
    ann[, .(cell_type, temporal_label, Notch, ntype)],
    by.x = by_x,
    by.y = "cell_type"
  )
  setorder(coord, .row_id)
  coord[, .row_id := NULL]
  
  return(coord)
}

#' Filter coordinates by gene expression data
filter_geneexp <- function(coord, ann, ts, syn_type = "pre") {
  by_x <- ifelse(syn_type == "pre", "pre_type", "post_type")

  # Extract gene symbols and convert to logical matrix
  ts_symbols <- ts$V1
  ts_data <- ts[, -1]  # Remove first column (gene names)

  # Annotate
  ann <- ann[Confident_annotation == "Y"]
  type_lut <- ann$cell_type[!is.na(ann$ozel2021_cluster)]
  names(type_lut) <- ann$ozel2021_cluster[!is.na(ann$ozel2021_cluster)]

  ts_data <- ts_data[, colnames(ts_data) %in% names(type_lut), with = FALSE]
  colnames(ts_data) <- type_lut[colnames(ts_data)]

  # Convert to logical and transpose to have cell types as rows
  ts_mat <- t(ts_data)
  colnames(ts_mat) <- ts_symbols

  # Keep only genes with at least one TRUE value
  ts_mat <- ts_mat[, colSums(ts_mat) > 0, drop = FALSE]
  ts_mat <- as.data.frame(ts_mat)
  ts_mat$cell_type <- row.names(ts_mat)

  # Annotate synapses with gene expression status
  coord[, .row_id := .I]
  coord <- merge(coord, ts_mat, by.x = by_x, by.y = "cell_type")

  coord <- merge(
    coord,
    ann[, .(cell_type, func, Notch, newly_ann, ntype)],
    by.x = by_x,
    by.y = "cell_type"
  )
  setorder(coord, .row_id)
  coord[, .row_id := NULL]

  return(coord)
}

#' Filter coordinates for LC/LPLC upstream partners
#' @param coord Synapse coordinate data (from rotated matrices)
#' @param ann Annotation data (visual_neurons_anno.csv)
#' @param partner_summary Presynaptic partner data (lc_lplc_presynaptic_partners.csv)
#' @param target_lc_lplc Specific LC/LPLC type to filter for
#' @param syn_type "pre" or "post" (default "pre" for upstream)
#' @return Filtered coordinates with contribution_score column
filter_lc_lplc_upstream <- function(coord, ann, partner_summary,
                                     target_lc_lplc, syn_type = "pre") {
  # Filter for target LC/LPLC type
  target_partners <- partner_summary[post_type == target_lc_lplc]

  # Calculate contribution scores (connectivity fraction)
  target_partners[, contribution_score := connected_neuron_count / total_neuron_count]

  # Create lookup: pre_type -> contribution_score
  score_lookup <- target_partners[, .(pre_type, contribution_score)]

  # Determine which column to filter on
  type_col <- paste0(syn_type, "_type")

  # Filter coordinates to keep only upstream partners
  coord <- coord[get(type_col) %in% score_lookup$pre_type]

  # Add row ID for order preservation
  coord[, .row_id := .I]

  # Merge with contribution scores
  coord <- merge(
    coord,
    score_lookup,
    by.x = type_col,
    by.y = "pre_type",
    all.x = TRUE
  )

  # Merge with annotations (for compatibility with existing viz system)
  if (nrow(ann) > 0) {
    coord <- merge(
      coord,
      ann[, .(cell_type, func, temporal_label, Notch, ntype)],
      by.x = type_col,
      by.y = "cell_type",
      all.x = TRUE
    )
  }

  # Restore order and clean up
  setorder(coord, .row_id)
  coord[, .row_id := NULL]

  return(coord)
}
