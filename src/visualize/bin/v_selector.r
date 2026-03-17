#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

options("ggrastr.default.dpi" = 450)

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

if (interactive()) {
  argvs$np <- "LOP_L"
  argvs$syn_type <- "post"
  argvs$ann <- "data/visual_neurons_anno.csv"
  argvs$ts <- "data/P15_CAM.csv"
  argvs$meta <- "data/viz_meta.csv"
  argvs$density <- "asis"
  argvs$subsample <- 10000L
  argvs$sparse_limit <- 100L
  source("./src/utils.r")
  syn_path <- file.path("int/idv_mat/", paste0(argvs$np, "_rotated.csv.gz"))
} else {
  argvs$subsample <- as.integer(argvs$subsample)
  argvs$sparse_limit <- as.integer(argvs$sparse_limit)
  syn_path <- argvs$synf
}


## Load plot metadata and presets with validation
plot_meta_all <- validate_config_file(
  argvs$meta,
  c(
    "neuropil", "x_axis", "y_axis", "axis_1_func", "axis_2_func", "min1",
    "max1", "min2", "max2", "zid", "zinvert", "minz", "maxz", "outlayout",
    "outr1", "outr2", "outd1", "outd2"
  )
)
plot_meta <- plot_meta_all[neuropil == argvs$np]
if (nrow(plot_meta) == 0) stop(sprintf("No metadata found for neuropil: %s", argvs$np))

# Generated 
x_axis <- paste(argvs$syn_type, plot_meta$x_axis, sep = "_")
y_axis <- paste(argvs$syn_type, plot_meta$y_axis, sep = "_")
scale_axis_1 <- get(plot_meta$axis_1_func)
scale_axis_2 <- get(plot_meta$axis_2_func)
color_func <- function() {
  f <- scale_color_manual(
    values = c(
      "Notch On" = "#008800", "Notch Off" = "#990099", "nil" = "grey75"
    ),
    breaks = c("Notch On", "Notch Off")
  )
  return(f)
}

## Load annotations, coordinates, and terminal selector expression
## (Ozel et al., 2022)
opc_anno <- fread(argvs$ann)
np_coord <- fread(syn_path)
ts <- fread(argvs$ts, header = TRUE)


## Validate that all data points fall within specified axis limits
x_values <- np_coord[[x_axis]]
y_values <- np_coord[[y_axis]]

# Annotate and filter to keep only annotated synapses
np_coord <- filter_type(
  np_coord, syn_type = argvs$syn_type, sparse_limit = argvs$sparse_limit
)
npp <- filter_geneexp(np_coord, opc_anno, ts, syn_type = argvs$syn_type)
np_plot <- rsubsample(npp, n = argvs$subsample, syn_type = argvs$syn_type)

ts_symbols <- colnames(np_plot)[
  # Assuming the only logical columns are all selector statuses
  vapply(np_plot, is.logical, logical(1))
]

# Generate shared legend once before the loop
temp_plot <- np_plot |>
  ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
  rasterize(geom_point(aes(color = Notch))) +
  labs(color = "Notch Status") +
  theme_ih2025() +
  color_func()

legendsp <- get_plot_component(temp_plot, "guide-box", return_all = TRUE)
# if (is.list(legendsp)) {
#   to_keep <- sapply(legendsp, function(x) "gtable" %in% class(x))
#   legendsp <- legendsp[to_keep]
#   if (length(legendsp) == 1) {
#     legendsp <- legendsp[[1]]
#   } else {
#     stop("Legend extraction error: There is more than 1 item.")
#   }
# }

for (i in ts_symbols) {
  if (!any(npp[[i]])) {
    next
  }
  
  out_prefix <- paste(
    argvs$np, argvs$syn_type, argvs$density, i,
    sep = "_"
  )
  
  # Create color labels for Notch On plot
  np_plot$on_label <- "nil"
  np_plot$on_label[np_plot[[i]] & np_plot$Notch == "Notch On"] <- "Notch On"
  
  # Create color labels for Notch Off plot
  np_plot$off_label <- "nil"
  np_plot$off_label[np_plot[[i]] & np_plot$Notch == "Notch Off"] <- "Notch Off"
  
  # Check if we have any non-nil labels to plot
  has_notch_on <- any(np_plot$on_label != "nil")
  has_notch_off <- any(np_plot$off_label != "nil")
  
  if (!has_notch_on && !has_notch_off) {
    message(sprintf("Skipping %s: no expressing synapses in either Notch population", i))
    next
  }
  
  if (grepl("reverse", plot_meta$axis_1_func)) {
    min1 <- plot_meta$max1
    max1 <- plot_meta$min1
    plot_meta$min1 <- min1
    plot_meta$max1 <- max1
  }
  if (grepl("reverse", plot_meta$axis_2_func)) {
    min2 <- plot_meta$max2
    max2 <- plot_meta$min2
    plot_meta$min2 <- min2
    plot_meta$max2 <- max2
  }
  
  
  # Generate plots only if they have data
  if (has_notch_on) {
    ## Generate Notch On plot
    dotp_on <- np_plot |>
      ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
      rasterize(geom_point(aes(color = on_label))) +
      labs(color = "Notch Status") +
      theme_ih2025() +
      color_func() +
      scale_axis_1(limits = c(plot_meta$min1, plot_meta$max1)) +
      scale_axis_2(limits = c(plot_meta$min2, plot_meta$max2)) +
      theme(legend.position="none")
  }
  
  if (has_notch_off) {
    ## Generate Notch Off plot
    dotp_off <- np_plot |>
      ggplot(aes(x = .data[[x_axis]], y = .data[[y_axis]])) +
      rasterize(geom_point(aes(color = off_label))) +
      labs(color = "Notch Status") +
      theme_ih2025() +
      color_func() +
      scale_axis_1(limits = c(plot_meta$min1, plot_meta$max1)) +
      scale_axis_2(limits = c(plot_meta$min2, plot_meta$max2)) +
      theme(legend.position="none")
  }
  
  ## Manual calculation and interpolation of density
  np_den <- npp[ , .(
    group = Notch,
    type = get(paste(argvs$syn_type, "type", sep = "_")),
    depth = get(paste(argvs$syn_type, "rz", sep = "_"))
  )]
  np_den$group[!npp[[i]]] <- "nil"
  
  if (plot_meta$zinvert) {
    np_den$depth <- np_den$depth * -1
  }
  np_den <- split(np_den, np_den$group, drop = TRUE)
  
  den_grid <- seq(plot_meta$minz, plot_meta$maxz, length.out = 1024)
  if (argvs$density == "asis") {
    np_den <- lapply(np_den, function(x) return(density(x$depth)))
    interpolated <- lapply(names(np_den), function(type) {
      d <- np_den[[type]]
      y_interp <- approx(d$x, d$y, xout = den_grid, rule = 2)$y
      return(data.frame(x = den_grid, y = y_interp, type = type))
    })
    inter <- do.call(rbind, interpolated)
  } else {
    type_avg <- lapply(names(np_den), function(tn) {
      x <- np_den[[tn]]
      id_type <- split(x, x$type, drop = TRUE)
      den_type <- lapply(id_type, function(idt) return(density(idt$depth)))
      if (length(den_type) > 0) {
        interpolated <- lapply(names(den_type), function(type) {
          d <- den_type[[type]]
          y_interp <- approx(d$x, d$y, xout = den_grid, rule = 2)$y
          return(data.frame(x = den_grid, y_interp = y_interp, type = type))
        })
        inter <- do.call(rbind, interpolated)
        out <- data.frame(
          x = den_grid,
          y = tapply(inter$y_interp, inter$x, mean),
          type = tn
        )
      } else {
        out <- data.frame(
          x = den_grid,
          y = 0,
          type = tn
        )
      }
      return(out)
    })
    inter <- do.call(rbind, type_avg)
  }
  
  denp <- inter |>
    ggplot(aes(x = x, y = y, color = type)) +
    geom_line() +
    theme_ih2025() +
    color_func() +
    guides(color = "none")
  
  if (plot_meta$zid == "y") denp <- denp + coord_flip()
  
  # Create and save combined plots only for populations with data
  if (has_notch_on) {
    # Create combined plots for Notch On
    if (plot_meta$outlayout == "landscape") {
      outp_on <- (dotp_on | denp) +
        plot_layout(widths = c(plot_meta$outr1, plot_meta$outr2))
    } else {
      outp_on <- (dotp_on / denp) +
        plot_layout(heights = c(plot_meta$outr1, plot_meta$outr2))
    }
    outp_on <- outp_on +
      plot_annotation(
        title = paste(i, "- Notch On"),
        subtitle = paste0(argvs$syn_type, "synapses")
      ) & 
      theme(
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "italic")
      )
    
    ggsave(
      plot = outp_on, filename = paste0(out_prefix, "_NotchOn.pdf"),
      width = plot_meta$outd1, height = plot_meta$outd2
    )
  }
  
  if (has_notch_off) {
    # Create combined plots for Notch Off
    if (plot_meta$outlayout == "landscape") {
      outp_off <- (dotp_off | denp) +
        plot_layout(widths = c(plot_meta$outr1, plot_meta$outr2))
    } else {
      outp_off <- (dotp_off / denp) +
        plot_layout(heights = c(plot_meta$outr1, plot_meta$outr2))
    }
    outp_off <- outp_off +
      plot_annotation(
        title = paste(i, "- Notch Off"),
        subtitle = paste0(argvs$syn_type, "synapses")
      ) & 
      theme(
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "italic")
      )
    
    ggsave(
      plot = outp_off, filename = paste0(out_prefix, "_NotchOff.pdf"),
      width = plot_meta$outd1, height = plot_meta$outd2
    )
  }
}

# Save shared legend
legend_prefix <- paste(
  argvs$np, argvs$syn_type, argvs$density,
  sep = "_"
)
ggsave(
  plot = legendsp, paste0(legend_prefix, "_legend.pdf")
)
