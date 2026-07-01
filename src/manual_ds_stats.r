library(ggplot2)
library(dplyr)
library(patchwork)
library(openxlsx)
library(RColorBrewer)
source("src/utils.r")

opc_anno <- read.csv("data/visual_neurons_anno.csv")
boot_sheets_paths <- "int/stats/deep_superficial/combined_broad_depth_results.xlsx"
# openxlsx::getSheetNames(boot_sheets_paths)
## [1] "broad_known"     "broad_new"       "Subsystem Known" "subsystem_new"   
##     "temporal_all"    "Temporal Known"  "temporal_new"    "type_putative"
boot <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "Temporal Known")
)

boot$types_of_interest <- factor(
  boot$types_of_interest,
  levels = c("Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", "Ey/Hbn", "Hbn/Opa/Slp", "Slp/D", "D/BH-1")
)

boot$syn_type <- factor(
  boot$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

ds_ribbon_alpha <- 0.1

# General plotting
dsplot <- function(
    stats, neuropil, split, title, bg_alpha = ds_ribbon_alpha,
    star_size = 4, star_location = 1.15,
    ylab = "Deep / Superficial Bias\n(+: Distal / -: Proximal)",
    scale_point_size_by_synapses = FALSE,
    synapse_size_range = c(0.8, 3.2),
    synapse_size_limits = NULL,
    synapse_size_breaks = waiver(),
    jitter_point_ranges = FALSE,
    point_jitter_width = 0.18,
    show_significance_stars = TRUE,
    synapse_point_stats = NULL
  ) {
  glossary <- c(
    "Notch Off_intrinsic" = "'Notch'^'Off'~'Interneurons'",
    "Notch Off_projection" = "'Notch'^'Off'~'Projection Neurons'",
    "Notch On_intrinsic" = "'Notch'^'On'~'Interneurons'",
    "Notch On_projection" = "'Notch'^'On'~'Projection Neurons'",
    "all" = "'Putative'~'OPC'~'Neurons'"
  )

  np_lut <- c(
    "ME_R" = "'Medulla'~",
    "LOP_R" = "'Lobula Plate'~",
    "LO_R" = "'Lobula'~"
  )
  spatial_groups <- c(
    "Vsx", "Optix", "Dpp",
    "Vsx/Optix", "Vsx/Dpp", "Optix/Dpp",
    "Vsx/Optix/Dpp"
  )
  uval <- unique(stats$types_of_interest[stats$neuropil == neuropil])
  is_spatial <- any(uval %in% spatial_groups)
  plot_title <- if (missing(title)) {
    paste0(
      np_lut[[neuropil]],
      ifelse(
        length(split) == 1,
        ifelse(is_spatial && split == "all", "'Spatial'~'Origin'~'Neurons'", glossary[[split]]),
        "'Projection'~'Neurons'"
      )
    )
  } else {
    title
  }
  broad_temporal_colors <- ih2025_colors()
  early_col <- broad_temporal_colors[["Early"]]
  late_col <- broad_temporal_colors[["Late"]]
  deep_col <- ifelse(grepl("^ME", neuropil), early_col, late_col)
  sup_col <- ifelse(grepl("^ME", neuropil), late_col, early_col)
  plot_data <- stats |>
    filter(neuropil == {{ neuropil }}, split %in% {{ split }})
  synapse_point_data <- NULL
  if (!is.null(synapse_point_stats)) {
    synapse_point_data <- synapse_point_stats |>
      filter(neuropil == {{ neuropil }}, split %in% {{ split }})
  }
  if (is_spatial) {
    plot_data$types_of_interest <- factor(
      plot_data$types_of_interest,
      levels = spatial_groups[spatial_groups %in% unique(plot_data$types_of_interest)]
    )
  }
  if (jitter_point_ranges) {
    set.seed(1)
    plot_data$.x_position <- as.numeric(plot_data$types_of_interest)
    if (!is.null(synapse_point_data)) {
      synapse_point_data$.x_position <- as.numeric(synapse_point_data$types_of_interest) +
        runif(nrow(synapse_point_data), -point_jitter_width, point_jitter_width)
    } else {
      plot_data$.x_position <- plot_data$.x_position +
        runif(nrow(plot_data), -point_jitter_width, point_jitter_width)
    }
  }

  pad_width <- max(
    length(unique(c(
      as.character(plot_data$types_of_interest),
      as.character(synapse_point_data$types_of_interest)
    ))) + 1,
    ifelse(is.factor(plot_data$types_of_interest), nlevels(plot_data$types_of_interest), 0) + 1
  )
  
  if (scale_point_size_by_synapses && is.null(synapse_point_data) && !"n_samples_interest" %in% names(plot_data)) {
    stop("Cannot scale point size: n_samples_interest column is missing")
  }
  if (scale_point_size_by_synapses && !is.null(synapse_point_data) && !"n_samples_interest" %in% names(synapse_point_data)) {
    stop("Cannot scale point size: n_samples_interest column is missing")
  }
  
  p <- plot_data %>%
    ggplot(aes(
      x = if (jitter_point_ranges) .data$.x_position else types_of_interest,
      color = types_of_interest
    )) +
    annotate(
      geom = "rect",
      fill = deep_col,
      ymin = -Inf,
      ymax = -0.5,
      xmin = 0,
      xmax = pad_width,
      alpha = bg_alpha
    ) +
    annotate(
      geom = "rect",
      fill = sup_col,
      ymin = 0.5,
      ymax = Inf,
      xmin = 0,
      xmax = pad_width,
      alpha = bg_alpha
    ) +
    geom_hline(yintercept = 0)

  if (scale_point_size_by_synapses && is.null(synapse_point_data)) {
    p <- p +
      geom_linerange(aes(
        ymin = bootstrap_bias_ratio_lower,
        ymax = bootstrap_bias_ratio_upper
      ), linewidth = 0.5) +
      geom_point(aes(
        y = bootstrap_bias_ratio_median,
        size = n_samples_interest
      )) +
      scale_size_continuous(
        trans = "log10",
        range = synapse_size_range,
        limits = synapse_size_limits,
        breaks = synapse_size_breaks,
        name = "Synapses",
        labels = scales::label_number()
      )
  } else {
    p <- p +
      geom_pointrange(aes(
        y = bootstrap_bias_ratio_median,
        ymin = bootstrap_bias_ratio_lower,
        ymax = bootstrap_bias_ratio_upper
      ), size = 0.5)
  }

  if (scale_point_size_by_synapses && !is.null(synapse_point_data)) {
    synapse_point_data <- synapse_point_data |>
      filter(
        !is.na(bootstrap_bias_ratio_median),
        !is.na(n_samples_interest),
        n_samples_interest > 0
      )
    p <- p +
      geom_point(
        data = synapse_point_data,
        aes(
          x = if (jitter_point_ranges) .data$.x_position else types_of_interest,
          y = bootstrap_bias_ratio_median,
          color = types_of_interest,
          size = n_samples_interest
        ),
        inherit.aes = FALSE,
        shape = 1,
        stroke = 0.6,
        alpha = 0.9
      ) +
      scale_size_continuous(
        trans = "log10",
        range = synapse_size_range,
        limits = synapse_size_limits,
        breaks = synapse_size_breaks,
        name = "Synapses",
        labels = scales::label_number()
      )
  }

  if (show_significance_stars) {
    p <- p + geom_text(
      aes(label = ifelse(.data$significant_fdr, "*", "")),
      y = star_location,
      color = "black",
      size = star_size
    )
  }

  p <- p +
    theme(
      panel.background = element_blank(),
      strip.background = element_rect(fill = "transparent", color = "black"),
      strip.text = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold")
    ) +
    guides(color = "none", size = guide_legend(order = 1)) +
    labs(
      y = ylab,
      title = parse(text = plot_title)
    ) +
    scale_y_continuous(limits = c(-1, star_location + 0.1))

  if (jitter_point_ranges) {
    type_levels <- levels(plot_data$types_of_interest)
    p <- p + scale_x_continuous(
      breaks = seq_along(type_levels),
      labels = type_levels,
      limits = c(0, pad_width),
      expand = expansion(mult = 0)
    )
  } else {
    p <- p + scale_x_discrete(drop = FALSE)
  }
  
  tws <- c("Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", "Ey/Hbn", "Hbn/Opa/Slp", "Slp/D", "D/BH-1")
  is_temporal <- any(uval %in% tws)
  
  sbs <- c("Color", "Form", "Luminance", "Motion", "Object", "Polarization", "Unannotated")
  is_subsystem <- any(uval %in% sbs)
  
  if (length(split) == 1) {
    if (split != "all") {
      if (is_temporal) {
        p <- p + scale_color_nk2023()
      } else {
        p <- p + scale_color_subsystem()
      }
    } else {
      if (is_subsystem) {
        p <- p + scale_color_subsystem()
      } else if (is_spatial) {
        p <- p + scale_color_spatial_origin()
      } else {
        p <- p + scale_color_type()
      }
      
    }
    return(p + facet_grid(~syn_type))
  } else {
    if (is_temporal) {
      p <- p + scale_color_nk2023()
    } else {
      p <- p + scale_color_subsystem()
    }
    return(p + facet_grid(split ~ syn_type))
  }
}
## Y1A — Medulla NotchOn interneurons: temporal shift
Y1A <- dsplot(
  boot, "ME_R", "Notch On_intrinsic",
)
Y1A  

## Y1B — Medulla NotchOff interneurons: temporal shift
Y1B <- dsplot(
  boot, "ME_R", "Notch Off_intrinsic",
)
Y1B  

# Y1C — Medulla NotchOn projection neurons: presyn temporal effect; postsyn distal
Y1C <- dsplot(
  boot, "ME_R", "Notch On_projection",
)
Y1C

# Y1D — Medulla NotchOff projection neurons
Y1D <- dsplot(
  boot, "ME_R", "Notch Off_projection",
)
Y1D

# Y1E — Lobula NotchOn projection neurons: early superficial → late deep
Y1E <- dsplot(
  boot, "LO_R", "Notch On_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1E

# Y1F — Lobula NotchOff projection neurons: no monotonic targeting
Y1F <- dsplot(
  boot, "LO_R", "Notch Off_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1F

# Y1G - Lobula plate TmY neurons: broadly distributed; exception TmY3 postsyn
Y1G <- dsplot(
  boot, "LOP_R", c("Notch Off_projection", "Notch On_projection"),
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1G

## New ones
boot_new <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "temporal_new")
)
boot_new$types_of_interest <- factor(
  boot_new$types_of_interest,
  levels = c("Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", "Ey/Hbn", "Hbn/Opa/Slp", "Slp/D", "D/BH-1")
)

boot_new$syn_type <- factor(
  boot_new$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

# Y1H - New Sm interneurons: weak/ambiguous bias
Y1H <- dsplot(
  boot_new, "ME_R", "Notch Off_intrinsic",
)
Y1H  

# Y1I - New lobula NotchOn Hbn/Opa/Slp types: deep-biased
Y1I <- dsplot(
  boot_new, "LO_R", "Notch On_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1I

# Y1J - Same cohort in medulla: subtly distal-biased
Y1J <- dsplot(
  boot_new, "ME_R", "Notch On_projection"
)
Y1J

# Y1K - New Erm/Ey NotchOff projection neurons in lobula: deep-biased
Y1K <- dsplot(
  boot_new, "LO_R", "Notch Off_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y1K

# Y1L - Same types in medulla: around serpentine; no prox/dist bias
Y1L <- dsplot(
  boot_new, "ME_R", "Notch Off_projection"
)
Y1L

Y1 <- list(Y1A, Y1B, Y1C, Y1D, Y1E, Y1F, Y1G, Y1H, Y1I, Y1J, Y1K, Y1L)

Y1_p <- wrap_plots(Y1, ncol = 3) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9))

ggsave(filename = "int/Supp_fig_Y1.pdf", plot = Y1_p, width = 8.5, height = 11)

type_depth <- read.csv("int/de_analysis/combined_type_depth.csv")
type_depth <- type_depth |>
  left_join(
    opc_anno |> select(cell_type, newly_ann),
    by = "cell_type"
  ) |>
  filter(
    !is.na(temporal_label),
    temporal_label != "",
    temporal_label != "unknown",
    is.na(skip_reason) | skip_reason == ""
  ) |>
  mutate(
    types_of_interest = factor(
      temporal_label,
      levels = c("Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", "Ey/Hbn", "Hbn/Opa/Slp", "Slp/D", "D/BH-1")
    ),
    split = paste(Notch, ntype, sep = "_"),
    n_samples_interest = n_syn_interest,
    syn_type = factor(
      syn_type,
      levels = c("pre", "post"),
      labels = c("Presynapse", "Postsynapse")
    )
  )

type_depth_known <- type_depth |>
  filter(is.na(newly_ann) | newly_ann != "Y")

type_depth_new <- type_depth |>
  filter(newly_ann == "Y")

y1_alt_synapse_size_limits <- range(type_depth$n_samples_interest, na.rm = TRUE)
y1_alt_synapse_size_breaks <- scales::breaks_log(n = 4)(y1_alt_synapse_size_limits)

Y1A_alt <- dsplot(
  boot, "ME_R", "Notch On_intrinsic",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_known
)
Y1B_alt <- dsplot(
  boot, "ME_R", "Notch Off_intrinsic",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_known
)
Y1C_alt <- dsplot(
  boot, "ME_R", "Notch On_projection",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_known
)
Y1D_alt <- dsplot(
  boot, "ME_R", "Notch Off_projection",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_known
)
Y1E_alt <- dsplot(
  boot, "LO_R", "Notch On_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_known
)
Y1F_alt <- dsplot(
  boot, "LO_R", "Notch Off_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_known
)
Y1G_alt <- dsplot(
  boot, "LOP_R", c("Notch Off_projection", "Notch On_projection"),
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_known
)
Y1H_alt <- dsplot(
  boot_new, "ME_R", "Notch Off_intrinsic",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_new
)
Y1I_alt <- dsplot(
  boot_new, "LO_R", "Notch On_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_new
)
Y1J_alt <- dsplot(
  boot_new, "ME_R", "Notch On_projection",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_new
)
Y1K_alt <- dsplot(
  boot_new, "LO_R", "Notch Off_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_new
)
Y1L_alt <- dsplot(
  boot_new, "ME_R", "Notch Off_projection",
  scale_point_size_by_synapses = TRUE,
  synapse_size_limits = y1_alt_synapse_size_limits,
  synapse_size_breaks = y1_alt_synapse_size_breaks,
  jitter_point_ranges = TRUE,
  synapse_point_stats = type_depth_new
)

Y1_alt <- list(
  Y1A_alt, Y1B_alt, Y1C_alt, Y1D_alt, Y1E_alt, Y1F_alt,
  Y1G_alt, Y1H_alt, Y1I_alt, Y1J_alt, Y1K_alt, Y1L_alt
)

Y1_alt_p <- wrap_plots(Y1_alt, ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 9),
    legend.position = "bottom"
  )

ggsave(filename = "int/Supp_fig_Y1-alt.pdf", plot = Y1_alt_p, width = 8.5, height = 11)

### Y2
putative_hl_cols <- grep("^putative_hl[0-9]+$", names(opc_anno), value = TRUE)
putative_hl_cols <- putative_hl_cols[
  colSums(opc_anno[putative_hl_cols] == "Y", na.rm = TRUE) > 0
]
putative_hl_ids <- as.integer(sub("^putative_hl", "", putative_hl_cols))
putative_hl_ids <- sort(putative_hl_ids)

boot_sheet_names <- getSheetNames(boot_sheets_paths)

read_putative_boot <- function(hl_id) {
  sheet_name <- paste0("type_putative_", hl_id)
  sheet_idx <- which(boot_sheet_names == sheet_name)
  if (length(sheet_idx) != 1) {
    stop(sprintf("Cannot find required worksheet: %s", sheet_name))
  }

  boot_putative <- read.xlsx(
    boot_sheets_paths,
    sheet = sheet_idx
  )

  boot_putative$syn_type <- factor(
    boot_putative$syn_type,
    levels = c("pre", "post"),
    labels = c("Presynapse", "Postsynapse")
  )

  boot_putative
}

make_y2_panel <- function(boot_putative, neuropil) {
  if (!any(boot_putative$neuropil == neuropil)) {
    return(NULL)
  }

  dsplot(
    boot_putative, neuropil, "all",
    ylab = ifelse(
      grepl("^LO", neuropil),
      "Deep Superficial Bias\n(+: Superficial / -: Deep)",
      "Deep / Superficial Bias\n(+: Distal / -: Proximal)"
    )
  )
}

primary_y2_panel_specs <- data.frame(
  hl_id = c(1, 2, 1, 2, 3),
  neuropil = c("LO_R", "LO_R", "ME_R", "ME_R", "ME_R"),
  stringsAsFactors = FALSE
)
primary_y2_panel_specs <- primary_y2_panel_specs[
  primary_y2_panel_specs$hl_id %in% putative_hl_ids,
]

additional_y2_panel_specs <- do.call(rbind, lapply(
  setdiff(putative_hl_ids, 1:3),
  function(hl_id) {
    data.frame(
      hl_id = hl_id,
      neuropil = c("LO_R", "ME_R"),
      stringsAsFactors = FALSE
    )
  }
))

y2_panel_specs <- rbind(primary_y2_panel_specs, additional_y2_panel_specs)

y2_boot <- lapply(
  putative_hl_ids,
  read_putative_boot
)
names(y2_boot) <- as.character(putative_hl_ids)

Y2 <- lapply(seq_len(nrow(y2_panel_specs)), function(i) {
  spec <- y2_panel_specs[i, ]
  make_y2_panel(y2_boot[[as.character(spec$hl_id)]], spec$neuropil)
})
Y2 <- Filter(Negate(is.null), Y2)

Y2_p <- wrap_plots(Y2, ncol = 2) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9))

ggsave(filename = "int/Supp_fig_Y2.pdf", plot = Y2_p, width = 12, height = 15)

### Y3
boot_fun <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "Subsystem Known")
)

boot_fun$syn_type <- factor(
  boot_fun$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

Y3A <- dsplot(
  boot_fun, "ME_R", "Notch On_intrinsic"
)
Y3A

Y3B <- dsplot(
  boot_fun, "ME_R", "Notch Off_intrinsic"
)
Y3B

Y3C <- dsplot(
  boot_fun, "LO_R", "Notch On_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y3C

Y3D <- dsplot(
  boot_fun, "LO_R", "Notch Off_projection",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y3D

Y3E <- dsplot(
  boot_fun, "ME_R", "Notch On_projection"
)
Y3E

Y3F <- dsplot(
  boot_fun, "ME_R", "Notch Off_projection"
)
Y3F

boot_fun_new <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "subsystem_new")
)

boot_fun_new$syn_type <- factor(
  boot_fun_new$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

Y3G <- dsplot(
  boot_fun_new, "ME_R", "Notch Off_intrinsic"
)
Y3G

Y3H <- dsplot(
  boot_fun_new, "ME_R", "Notch Off_projection"
)
Y3H

Y3I <- dsplot(
  boot_fun_new, "ME_R", "Notch On_projection"
)
Y3I

Y3J <- dsplot(
  boot_fun_new, "LO_R", c("Notch Off_projection", "Notch On_projection"),
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y3J

Y3K <- dsplot(
  boot_fun, "LOP_R", c("Notch On_projection", "Notch Off_projection"),
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)

Y3K

Y3 <- list(Y3A, Y3B, Y3C, Y3D, Y3E, Y3F, Y3G, Y3H, Y3I, Y3J, Y3K)

Y3_p <- wrap_plots(Y3, ncol = 3) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9))

ggsave(filename = "int/Supp_fig_Y3.pdf", plot = Y3_p, width = 8.5, height = 11)

### Y4
boot_fun_putative <- read.xlsx(
  boot_sheets_paths,
  sheet = which(getSheetNames(boot_sheets_paths) == "subsystem_putative")
)

boot_fun_putative$syn_type <- factor(
  boot_fun_putative$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

Y4A <- dsplot(
  boot_fun_putative, "ME_R", "all"
)
Y4A

Y4B <- dsplot(
  boot_fun_putative, "LO_R", "all",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y4B

Y4C <-dsplot(
  boot_fun_putative, "LOP_R", "all",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
) 
Y4C

Y4 <- list(
  Y4A, Y4B, Y4C,
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer(),
  plot_spacer()
)

Y4_p <- wrap_plots(Y4, ncol = 3) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9))

ggsave(filename = "int/Supp_fig_Y4.pdf", plot = Y4_p, width = 8.5, height = 11)

### Y6
spatial_sheet_idx <- which(getSheetNames(boot_sheets_paths) == "spatial_all")
if (length(spatial_sheet_idx) != 1) {
  stop("Cannot find required worksheet: spatial_all")
}

boot_spatial <- read.xlsx(
  boot_sheets_paths,
  sheet = spatial_sheet_idx
)

boot_spatial$syn_type <- factor(
  boot_spatial$syn_type,
  levels = c("pre", "post"),
  labels = c("Presynapse", "Postsynapse")
)

Y6A <- dsplot(
  boot_spatial, "ME_R", "all"
)
Y6A

Y6B <- dsplot(
  boot_spatial, "LO_R", "all",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y6B

Y6C <- dsplot(
  boot_spatial, "LOP_R", "all",
  ylab = "Deep Superficial Bias\n(+: Superficial / -: Deep)"
)
Y6C

Y6 <- list(
  Y6A, Y6B, Y6C,
  plot_spacer(),
  plot_spacer(),
  plot_spacer()
)

Y6_p <- wrap_plots(Y6, ncol = 3) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 9))

ggsave(filename = "int/Supp_fig_Y6.pdf", plot = Y6_p, width = 8.5, height = 4)
