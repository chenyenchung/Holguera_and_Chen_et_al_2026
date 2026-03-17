#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Source utility functions (soft-linked by Nextflow)
if (file.exists("./utils.r")) {
  source("./utils.r", chdir = FALSE)
}

### TODO: Interactive mode for development
if (interactive()) {
  source("./src/utils.r", chdir = FALSE)
  argvs$partner_data <- "../Flywire_Isabel/partner_data/partner_summary.csv"
  argvs$threshold <- 0.03
  argvs$syn <- "both"
  argvs$np <- "All"
} else {
  argvs$threshold <- as.numeric(argvs$threshold)
}

# Hard-coded types of interest
types_of_interest <- c("Tm5f", "Tm8a", "Li06", "TmY10", "TmY11", "Sm03", "Tm5d",
                      "Tm36", "Tm37", "Sm01", "Sm02", "Sm14", "Tm5e", "Tm33")

# Load partner data
cat("Loading partner data...\n")
partner_data <- fread(argvs$partner_data)

# Function to analyze partners for a neuropil using synapse count
analyze_neuropil_partners_synapses <- function(
    data, neuropil_name, types_of_interest, syn_type, threshold
  ) {
  if (syn_type != "both") {
    syn_use <- paste0(syn_type, "_type")
    partner_use <- ifelse(syn_use == "pre_type", "post_type", "pre_type")
    if (neuropil_name != "All") {
      # Filter for types of interest as presynaptic partners in this neuropil
      connections <- data[
        neuropil == neuropil_name & get(syn_use) %in% types_of_interest
      ]
    } else {
      connections <- data[
        get(syn_use) %in% types_of_interest
      ]
    }
    
    setnames(connections, syn_use, "type")
    setnames(connections, partner_use, "partner")
  } else {
    connections_com <- copy(data)
    setnames(connections_com, "pre_type", "type")
    setnames(connections_com, "post_type", "partner")
    connections <- copy(data)
    setnames(connections, "post_type", "type")
    setnames(connections, "pre_type", "partner")
    connections <- rbind(connections, connections_com)
    
    if (neuropil_name != "All") {
      connections <- connections[
        neuropil == neuropil_name & type %in% types_of_interest
      ]
    } else {
      connections <- connections[type %in% types_of_interest]
    }
  }
  
  # Calculate ratios and filter by threshold
  partner_ratios <- connections[
    , ratio := nSynapse / sum(nSynapse),
    by = "type"
  ][, category := fifelse(ratio > argvs$threshold, partner, "Others")][
    , .(ratio = sum(ratio)), by = c("type", "category")
  ]
  
  return(partner_ratios)
}

# Function to create donut plot
create_donut_plot <- function(
    plot_data, neuropil_name, syn_type, threshold
    ) {
  # Pool all data across types of interest for the neuropil
  type_use <- paste0(syn_type, "_type")
  pooled_data <- plot_data
  setorder(pooled_data, type, ratio)
  
  this_type <- unique(pooled_data$type)
  if (length(this_type) > 1) {
    stop("The plot_data should only contain one type")
  }
  
  pooled_data[ , ymax := cumsum(ratio)][
    , ymin := c(0, head(ymax, n = -1))
  ][
    , y := (ymin + ymax) / 2
  ]
  
  # Create plot
  p <- ggplot(pooled_data, aes(
    ymin = ymin, ymax = ymax, group = category, 
    fill = ratio, xmin = 3, xmax = 12
  )) +
    geom_rect(color = "black") +
    geom_label_repel(
      aes(x = 9, y = y,
          label = paste0(category, ": ", round(ratio * 100, 1), "%")),
      nudge_x = 5, fill = "white", size = 2
    ) +
    scale_x_continuous(limits = c(0, 18)) +
    scale_fill_gradient(low = ifelse(syn_type == "pre", "lightblue", "pink"),
                        high = ifelse(syn_type == "pre", "blue", "#800020"),
                        limits = c(0, NA), labels = scales::percent) +
    labs(
      title = paste0("Main ", syn_type,"synaptic partners in ", neuropil_name),
      subtitle = paste0(this_type,
                        "\n(Partners accounting for < ",
                        threshold * 100,
                        "% are considered others)"),
      fill = "% of connections"
    ) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  
  return(p)
}

# Process each neuropil
cat(sprintf("Processing %s neuropil...\n", argvs$np))

# Analyze partners using synapse count
plot_data <- analyze_neuropil_partners_synapses(
  partner_data, argvs$np, types_of_interest, argvs$syn, argvs$threshold
)

plot_data <- split(plot_data, plot_data$type, drop = TRUE)

for (ntype in names(plot_data)) {
  # Create and save plot
  p <- create_donut_plot(
    plot_data[[ntype]], argvs$np, argvs$syn, argvs$threshold
  )
  
  output_file <- paste(
    argvs$np, ntype, argvs$syn, "partners.pdf", sep = "_"
  )
  ggsave(output_file, plot = p, width = 8, height = 6)
  cat(sprintf("Saved %s\n", output_file))
}

cat("Analysis complete!\n")