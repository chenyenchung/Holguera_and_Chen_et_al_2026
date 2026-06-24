#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(scales))

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }

    name <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[name]] <- TRUE
      i <- i + 1L
    } else {
      out[[name]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

argvs <- parse_args(commandArgs(trailingOnly = TRUE))
argvs$ann <- if (is.null(argvs$ann)) "data/visual_neurons_anno.csv" else argvs$ann
argvs$rotated_dir <- if (is.null(argvs$rotated_dir)) "int/idv_mat" else argvs$rotated_dir
argvs$output_dir <- if (is.null(argvs$output_dir)) "int/origin_heatmap" else argvs$output_dir
argvs$output_prefix <- if (is.null(argvs$output_prefix)) "origin_heatmap" else argvs$output_prefix

temporal_levels <- c(
  "Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey",
  "Ey/Hbn", "Hbn/Opa/Slp", "Slp/D", "D/BH-1"
)

spatial_families <- list(
  Vsx = c("vVsx", "dVsx"),
  Optix = c("vOptix", "dOptix"),
  Dpp = c("vDpp", "dDpp")
)
spatial_levels <- c(
  "Vsx", "Optix", "Dpp",
  "Vsx/Optix", "Vsx/Dpp", "Optix/Dpp", "Vsx/Optix/Dpp"
)

neuropils <- c("ME_R", "LO_R", "LOP_R")
syn_types <- c("pre", "post")

is_origin_flag <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  !is.na(x_chr) & x_chr %in% c("1", "TRUE", "T", "Y", "YES")
}

prepare_annotation <- function(path) {
  required_cols <- c(
    "cell_type", "temporal_label", "Confident_annotation",
    unlist(spatial_families, use.names = FALSE)
  )
  ann <- fread(path)
  missing_cols <- setdiff(required_cols, names(ann))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Missing required annotation columns in %s: %s",
      path,
      paste(missing_cols, collapse = ", ")
    ))
  }

  for (col in unlist(spatial_families, use.names = FALSE)) {
    ann[, (col) := is_origin_flag(get(col))]
  }

  family_order <- names(spatial_families)
  family_flag_cols <- paste0(".spatial_", family_order)
  names(family_flag_cols) <- family_order
  for (family in family_order) {
    ann[
      ,
      (family_flag_cols[[family]]) := rowSums(.SD) > 0,
      .SDcols = spatial_families[[family]]
    ]
  }

  ann[, spatial_origin := apply(.SD, 1, function(row) {
    origins <- family_order[as.logical(row)]
    if (length(origins) == 0) {
      return(NA_character_)
    }
    paste(origins, collapse = "/")
  }), .SDcols = family_flag_cols]

  ann <- ann[
    Confident_annotation == "Y" &
      temporal_label %in% temporal_levels &
      spatial_origin %in% spatial_levels,
    .(cell_type, temporal_origin = temporal_label, spatial_origin)
  ]

  ann[, temporal_origin := factor(temporal_origin, levels = temporal_levels)]
  ann[, spatial_origin := factor(spatial_origin, levels = spatial_levels)]
  ann
}

complete_origin_grid <- function() {
  CJ(
    temporal_origin = temporal_levels,
    spatial_origin = spatial_levels,
    unique = TRUE
  )
}

count_panel <- function(np, syn_type, ann, rotated_dir) {
  syn_path <- file.path(rotated_dir, paste0(np, "_rotated.csv.gz"))
  if (!file.exists(syn_path)) {
    stop(sprintf("Cannot find rotated synapse matrix: %s", syn_path))
  }

  type_col <- paste0(syn_type, "_type")
  coord <- fread(syn_path, select = type_col)
  if (!type_col %in% names(coord)) {
    stop(sprintf("Missing required column %s in %s", type_col, syn_path))
  }

  setnames(coord, type_col, "cell_type")
  counts <- merge(coord, ann, by = "cell_type", allow.cartesian = TRUE)[
    ,
    .(count = .N),
    by = .(temporal_origin, spatial_origin)
  ]

  grid <- complete_origin_grid()
  out <- merge(grid, counts, by = c("temporal_origin", "spatial_origin"), all.x = TRUE)
  out[is.na(count), count := 0L]
  out[, `:=`(
    neuropil = np,
    syn_type = syn_type,
    count_log10 = fifelse(count > 0, log10(count), NA_real_)
  )]
  setcolorder(out, c("neuropil", "syn_type", "temporal_origin", "spatial_origin", "count", "count_log10"))
  out[]
}

plot_panel <- function(panel_counts) {
  np <- unique(panel_counts$neuropil)
  syn_type <- unique(panel_counts$syn_type)
  if (length(np) != 1L || length(syn_type) != 1L) {
    stop("plot_panel expects counts from exactly one neuropil and synapse type")
  }

  panel_counts <- copy(panel_counts)
  panel_counts[, temporal_origin := factor(temporal_origin, levels = temporal_levels)]
  panel_counts[, spatial_origin := factor(spatial_origin, levels = rev(spatial_levels))]
  panel_counts[, fill_count := fifelse(count > 0, as.numeric(count), NA_real_)]

  max_count <- max(panel_counts$count, na.rm = TRUE)
  fill_limits <- c(1, max(1, max_count))

  ggplot(panel_counts, aes(x = temporal_origin, y = spatial_origin, fill = fill_count)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = ifelse(count > 0, count, "")), size = 2.4, color = "black") +
    scale_fill_gradient(
      low = "#F7FBFF",
      high = "#08519C",
      trans = "log10",
      limits = fill_limits,
      oob = squish,
      na.value = "grey95",
      labels = label_number(),
      name = "Synapse\ncount"
    ) +
    labs(
      title = paste(np, toupper(syn_type), sep = " "),
      x = "Temporal origin",
      y = "Spatial origin"
    ) +
    coord_equal() +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.key.height = grid::unit(0.35, "in")
    )
}

if (!file.exists(argvs$ann)) {
  stop(sprintf("Cannot find annotation file: %s", argvs$ann))
}
if (!dir.exists(argvs$rotated_dir)) {
  stop(sprintf("Cannot find rotated matrix directory: %s", argvs$rotated_dir))
}
if (!dir.exists(argvs$output_dir)) {
  dir.create(argvs$output_dir, recursive = TRUE)
}

anno <- prepare_annotation(argvs$ann)
if (nrow(anno) == 0) {
  stop("No confidently annotated origin rows remain after filtering")
}

counts <- rbindlist(lapply(neuropils, function(np) {
  rbindlist(lapply(syn_types, function(syn_type) {
    count_panel(np, syn_type, anno, argvs$rotated_dir)
  }))
}))

counts_path <- file.path(argvs$output_dir, paste0(argvs$output_prefix, "_counts.csv"))
fwrite(counts, counts_path)

plots <- lapply(neuropils, function(np) {
  lapply(syn_types, function(syn_type) {
    panel_syn_type <- syn_type
    plot_panel(counts[neuropil == np & syn_type == panel_syn_type])
  })
})
plots <- unlist(plots, recursive = FALSE)

combined <- wrap_plots(plots, ncol = 2)
pdf_path <- file.path(argvs$output_dir, paste0(argvs$output_prefix, ".pdf"))
ggsave(pdf_path, combined, width = 11, height = 10, units = "in")

message(sprintf("Wrote counts: %s", counts_path))
message(sprintf("Wrote heatmap: %s", pdf_path))
