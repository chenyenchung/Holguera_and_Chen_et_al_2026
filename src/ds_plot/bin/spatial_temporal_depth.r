#!/usr/bin/env Rscript
script_file <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
repo_root <- normalizePath(
  file.path(dirname(script_file), "../../.."),
  mustWork = TRUE
)
renv::load(repo_root)
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)
required_args <- c("ann", "utils", "type_depth", "out_dir")
missing_args <- required_args[vapply(
  required_args,
  function(x) is.null(argvs[[x]]),
  logical(1)
)]
if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

source(argvs$utils, chdir = FALSE)
dir.create(argvs$out_dir, recursive = TRUE, showWarnings = FALSE)

temporal_levels <- c(
  "Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey",
  "Ey/Hbn", "Hbn/Opa/Slp", "Slp/D", "D/BH-1"
)
spatial_levels <- c(
  "Vsx", "Optix", "Dpp", "Vsx/Optix", "Vsx/Dpp",
  "Optix/Dpp", "Vsx/Optix/Dpp"
)
stratum_levels <- c(
  "Notch On intrinsic",
  "Notch Off intrinsic",
  "Notch On projection",
  "Notch Off projection"
)
neuropils <- c("ME_R", "LO_R", "LOP_R")
panel_letters <- setNames(c("a", "b", "c"), neuropils)

required_annotation_cols <- c(
  "cell_type", "temporal_label", "temporal_id", "Notch", "ntype",
  "Confident_annotation", "newly_ann", "ozel2021_cluster",
  "vVsx", "dVsx", "vOptix", "dOptix", "vDpp", "dDpp"
)
annotation <- fread(argvs$ann)
missing_annotation_cols <- setdiff(required_annotation_cols, names(annotation))
if (length(missing_annotation_cols) > 0) {
  stop(
    "Missing annotation columns: ",
    paste(missing_annotation_cols, collapse = ", ")
  )
}

annotation <- add_spatial_origin(annotation)
spatial_table <- annotation[
  Confident_annotation == "Y" &
    temporal_label %in% temporal_levels &
    as.character(spatial_origin) %in% spatial_levels,
  .(
    cell_type,
    temporal_origin = temporal_label,
    temporal_id,
    spatial_origin = as.character(spatial_origin),
    notch_status = Notch,
    neuron_class = ntype,
    confident_annotation = Confident_annotation,
    newly_annotated = newly_ann,
    ozel2021_cluster,
    vVsx, dVsx, vOptix, dOptix, vDpp, dDpp
  )
]
if (nrow(spatial_table) == 0) {
  stop("No confidently annotated temporal/spatial neuronal types were found")
}
if (anyDuplicated(spatial_table$cell_type)) {
  stop("Spatial-origin table contains duplicated cell types")
}
spatial_table[, temporal_origin := factor(
  temporal_origin,
  levels = temporal_levels
)]
spatial_table[, spatial_origin := factor(
  spatial_origin,
  levels = spatial_levels
)]
setorder(spatial_table, temporal_origin, spatial_origin, cell_type)

required_depth_cols <- c(
  "cell_type", "Notch", "ntype", "temporal_label", "neuropil",
  "syn_type", "skip_reason", "direction", "significant_fdr", "p_value_fdr",
  "bootstrap_bias_ratio_median", "bootstrap_bias_ratio_lower",
  "bootstrap_bias_ratio_upper", "n_neurons_interest", "n_syn_interest"
)
depth <- fread(argvs$type_depth)
missing_depth_cols <- setdiff(required_depth_cols, names(depth))
if (length(missing_depth_cols) > 0) {
  stop("Missing type-depth columns: ", paste(missing_depth_cols, collapse = ", "))
}
depth[, skip_reason := trimws(as.character(skip_reason))]
depth[skip_reason == "", skip_reason := NA_character_]

figure_data <- merge(
  depth,
  spatial_table[, .(
    cell_type,
    mapped_temporal_origin = as.character(temporal_origin),
    spatial_origin = as.character(spatial_origin),
    mapped_notch_status = notch_status,
    mapped_neuron_class = neuron_class
  )],
  by = "cell_type",
  all = FALSE
)

metadata_mismatch <- figure_data[
  temporal_label != mapped_temporal_origin |
    Notch != mapped_notch_status |
    ntype != mapped_neuron_class
]
if (nrow(metadata_mismatch) > 0) {
  stop(
    "Type-depth metadata disagrees with the spatial annotation for: ",
    paste(unique(metadata_mismatch$cell_type), collapse = ", ")
  )
}

figure_data[, `:=`(
  temporal_origin = factor(mapped_temporal_origin, levels = temporal_levels),
  spatial_origin = factor(spatial_origin, levels = spatial_levels),
  notch_status = mapped_notch_status,
  neuron_class = mapped_neuron_class,
  developmental_stratum = factor(
    paste(mapped_notch_status, mapped_neuron_class),
    levels = stratum_levels
  ),
  syn_type = factor(
    syn_type,
    levels = c("pre", "post"),
    labels = c("Presynapse", "Postsynapse")
  )
)]
figure_data[, included_in_figure :=
  neuropil %in% neuropils &
    is.na(skip_reason) &
    !is.na(bootstrap_bias_ratio_median) &
    !is.na(developmental_stratum)]

# Give each neuronal type a stable horizontal offset. Spatial-origin offsets
# separate colors within a temporal window; the smaller type offset prevents
# exact overlap without making plot positions depend on row order.
stable_type_offset <- function(cell_type) {
  bytes <- utf8ToInt(enc2utf8(cell_type))
  hash <- sum(bytes * seq_along(bytes)) %% 101
  (hash / 100 - 0.5) * 0.05
}
figure_data[, x_position :=
  as.numeric(temporal_origin) +
    (as.numeric(spatial_origin) - 4) * 0.085 +
    vapply(cell_type, stable_type_offset, numeric(1))]

plot_data <- figure_data[included_in_figure == TRUE]
if (nrow(plot_data) == 0) {
  stop("No successful type-depth rows are available for the spatial figure")
}

stratum_title_lut <- c(
  "Notch On intrinsic" = "'Notch'^'On'~'Interneurons'",
  "Notch Off intrinsic" = "'Notch'^'Off'~'Interneurons'",
  "Notch On projection" = "'Notch'^'On'~'Projection Neurons'",
  "Notch Off projection" = "'Notch'^'Off'~'Projection Neurons'"
)
neuropil_title_lut <- c(
  ME_R = "'Medulla'~",
  LO_R = "'Lobula'~",
  LOP_R = "'Lobula Plate'~"
)
broad_temporal_colors <- ih2025_colors()
early_col <- broad_temporal_colors[["Early"]]
late_col <- broad_temporal_colors[["Late"]]
pad_width <- length(temporal_levels) + 1
neuropil_column_width <- 4
stratum_row_height <- 2.75

make_panel <- function(np, stratum) {
  current <- copy(plot_data[
    neuropil == np & as.character(developmental_stratum) == stratum
  ])
  if (nrow(current) == 0) {
    return(plot_spacer())
  }

  y_label <- if (np == "ME_R") {
    "Deep / Superficial Bias\n(+: Distal / -: Proximal)"
  } else {
    "Deep Superficial Bias\n(+: Superficial / -: Deep)"
  }
  deep_col <- if (np == "ME_R") early_col else late_col
  superficial_col <- if (np == "ME_R") late_col else early_col
  plot_title <- paste0(
    neuropil_title_lut[[np]],
    stratum_title_lut[[stratum]]
  )

  ggplot(current, aes(x = x_position, color = spatial_origin)) +
    annotate(
      "rect", xmin = 0, xmax = pad_width, ymin = -Inf, ymax = -0.5,
      fill = deep_col, alpha = 0.1
    ) +
    annotate(
      "rect", xmin = 0, xmax = pad_width, ymin = 0.5, ymax = Inf,
      fill = superficial_col, alpha = 0.1
    ) +
    geom_hline(yintercept = 0) +
    geom_pointrange(
      aes(
        y = bootstrap_bias_ratio_median,
        ymin = bootstrap_bias_ratio_lower,
        ymax = bootstrap_bias_ratio_upper
      ),
      size = 0.5
    ) +
    facet_grid(~syn_type) +
    scale_x_continuous(
      breaks = seq_along(temporal_levels),
      labels = temporal_levels,
      limits = c(0, pad_width),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_color_spatial_origin() +
    labs(
      title = parse(text = plot_title),
      x = NULL,
      y = y_label,
      color = "Spatial origin"
    ) +
    theme(
      panel.background = element_blank(),
      strip.background = element_rect(fill = "transparent", color = "black"),
      strip.text = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold"),
      plot.tag = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
}

panels <- setNames(lapply(neuropils, function(np) {
  setNames(lapply(stratum_levels, function(stratum) {
    make_panel(np, stratum)
  }), stratum_levels)
}), neuropils)

observed_spatial_levels <- spatial_levels[
  spatial_levels %in% unique(as.character(plot_data$spatial_origin))
]
legend_data <- data.table(
  spatial_origin = factor(
    observed_spatial_levels,
    levels = observed_spatial_levels
  ),
  x = seq_along(observed_spatial_levels),
  y = 1
)
legend_source <- ggplot(
  legend_data,
  aes(x = x, y = y, color = spatial_origin)
) +
  geom_point(size = 2) +
  scale_color_spatial_origin() +
  labs(color = "Spatial origin") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6)
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
spatial_origin_legend <- cowplot::get_legend(legend_source)
legend_path <- file.path(
  argvs$out_dir,
  "Supp_Figure_18_spatial_within_temporal_legend.pdf"
)
ggsave(
  legend_path,
  cowplot::ggdraw(spatial_origin_legend),
  width = 5.5,
  height = 0.65,
  units = "in"
)

make_neuropil_column <- function(np, tagged = FALSE, drop_empty = FALSE) {
  panel_strata <- stratum_levels
  if (drop_empty) {
    panel_strata <- panel_strata[vapply(
      panel_strata,
      function(stratum) {
        nrow(plot_data[
          neuropil == np & as.character(developmental_stratum) == stratum
        ]) > 0
      },
      logical(1)
    )]
  }
  current <- panels[[np]][panel_strata]
  if (tagged) {
    current[[1]] <- current[[1]] + labs(tag = panel_letters[[np]])
  }
  wrap_plots(current, ncol = 1) & theme(legend.position = "none")
}

for (np in neuropils) {
  panel_path <- file.path(
    argvs$out_dir,
    sprintf(
      "Supp_Figure_18%s_%s_spatial_within_temporal.pdf",
      panel_letters[[np]],
      np
    )
  )
  ggsave(
    panel_path,
    make_neuropil_column(np, drop_empty = TRUE),
    width = neuropil_column_width,
    height = stratum_row_height * length(unique(
      plot_data[neuropil == np, developmental_stratum]
    )),
    units = "in"
  )
}

review_plot <- wrap_plots(
  lapply(neuropils, function(np) make_neuropil_column(np, tagged = TRUE)),
  ncol = 3
) & theme(legend.position = "none")
review_path <- file.path(
  argvs$out_dir,
  "Supp_Figure_18_spatial_within_temporal_review.pdf"
)
ggsave(
  review_path,
  review_plot,
  width = neuropil_column_width * length(neuropils),
  height = stratum_row_height * length(stratum_levels),
  units = "in"
)

spatial_csv <- copy(spatial_table)
spatial_csv[, `:=`(
  temporal_origin = as.character(temporal_origin),
  spatial_origin = as.character(spatial_origin)
)]
spatial_csv_path <- file.path(
  argvs$out_dir,
  "spatial_origin_neuronal_types.csv"
)
fwrite(spatial_csv, spatial_csv_path)

figure_export <- figure_data[, .(
  cell_type,
  temporal_origin = as.character(temporal_origin),
  spatial_origin = as.character(spatial_origin),
  notch_status,
  neuron_class,
  neuropil,
  synapse_type = as.character(syn_type),
  included_in_figure,
  skip_reason,
  direction,
  significant_fdr,
  p_value_fdr,
  bootstrap_bias_ratio_median,
  bootstrap_bias_ratio_lower,
  bootstrap_bias_ratio_upper,
  n_neurons_interest,
  n_syn_interest
)]
setorder(
  figure_export,
  neuropil,
  synapse_type,
  temporal_origin,
  spatial_origin,
  cell_type
)

workbook_path <- file.path(
  argvs$out_dir,
  "spatial_origin_neuronal_types.xlsx"
)
wb <- createWorkbook()
addWorksheet(wb, "neuronal_types")
writeData(wb, "neuronal_types", spatial_csv, withFilter = TRUE)
freezePane(wb, "neuronal_types", firstRow = TRUE)
setColWidths(wb, "neuronal_types", cols = seq_len(ncol(spatial_csv)), widths = "auto")
addWorksheet(wb, "figure_data")
writeData(wb, "figure_data", figure_export, withFilter = TRUE)
freezePane(wb, "figure_data", firstRow = TRUE)
setColWidths(wb, "figure_data", cols = seq_len(ncol(figure_export)), widths = "auto")
saveWorkbook(wb, workbook_path, overwrite = TRUE)

message("Wrote spatial/temporal figure panels and review sheet to: ", argvs$out_dir)
message("Wrote spatial-origin CSV: ", spatial_csv_path)
message("Wrote spatial-origin workbook: ", workbook_path)
