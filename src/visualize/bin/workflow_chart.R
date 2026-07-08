#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(R.utils))

argvs <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

output <- argvs$output
if (is.null(output) || is.na(output) || !nzchar(output)) {
  output <- "int/workflow_chart/flyem_workflow_chart.pdf"
}

out_dir <- dirname(output)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

colors <- list(
  ink = "#1F2933",
  muted = "#5B6470",
  line = "#6B7280",
  lane = "#EEF2F6",
  source = "#E8F2FA",
  prep = "#E8F5ED",
  analysis = "#FFF4D8",
  output = "#F2EAF7",
  source_border = "#4C8CB5",
  prep_border = "#4C9A67",
  analysis_border = "#C9952D",
  output_border = "#8563A8"
)

nodes <- data.frame(
  id = c(
    "connectome", "annotations", "distances",
    "filter", "align",
    "maps", "depth", "function", "genes", "context",
    "figures", "tables"
  ),
  lane = c(
    "Data sources", "Data sources", "Data sources",
    "Preprocessing", "Preprocessing",
    "Analysis branches", "Analysis branches", "Analysis branches",
    "Analysis branches", "Analysis branches",
    "Outputs", "Outputs"
  ),
  title = c(
    "FlyWire visual-system connectome",
    "Developmental annotations and gene-expression resources",
    "Published connectivity-distance data",
    "Filter and annotate visual-system synapses",
    "Align neuropil coordinates",
    "Spatial synapse maps",
    "Superficial/deep depth bias",
    "Functional enrichment",
    "Gene-expression targeting",
    "Connectivity context",
    "Publication figures",
    "Statistical summaries and candidate tables"
  ),
  summary = c(
    "",
    "",
    "",
    "",
    "Split by neuropil and rotate synapse positions into a shared depth frame.",
    "Compare synapse locations by temporal origin, function, and cell type.",
    "Quantify whether groups preferentially target reference-defined depth domains.",
    "Test whether temporal and Notch-defined cohorts associate with functional subsystems.",
    "Relate selector and CAM expression programs to synapse depth and spatial maps.",
    "Summarize cell-type similarity and major synaptic partner structure.",
    "",
    ""
  ),
  x = c(
    0.135, 0.135, 0.135,
    0.365, 0.365,
    0.620, 0.620, 0.620, 0.620, 0.620,
    0.885, 0.885
  ),
  y = c(
    0.710, 0.505, 0.300,
    0.610, 0.390,
    0.790, 0.625, 0.460, 0.295, 0.130,
    0.610, 0.390
  ),
  w = c(
    0.205, 0.205, 0.205,
    0.205, 0.205,
    0.240, 0.240, 0.240, 0.240, 0.240,
    0.205, 0.205
  ),
  h = c(
    0.125, 0.145, 0.125,
    0.140, 0.155,
    0.120, 0.135, 0.145, 0.145, 0.135,
    0.130, 0.145
  ),
  fill = c(
    colors$source, colors$source, colors$source,
    colors$prep, colors$prep,
    colors$analysis, colors$analysis, colors$analysis, colors$analysis, colors$analysis,
    colors$output, colors$output
  ),
  border = c(
    colors$source_border, colors$source_border, colors$source_border,
    colors$prep_border, colors$prep_border,
    colors$analysis_border, colors$analysis_border, colors$analysis_border,
    colors$analysis_border, colors$analysis_border,
    colors$output_border, colors$output_border
  ),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c(
    "connectome", "annotations", "filter", "align", "align", "align", "align",
    "annotations", "annotations", "annotations", "distances",
    "maps", "depth", "function", "genes", "context"
  ),
  to = c(
    "filter", "filter", "align", "maps", "depth", "genes", "context",
    "maps", "depth", "function", "context",
    "figures", "figures", "tables", "figures", "figures"
  ),
  stringsAsFactors = FALSE
)

lane_labels <- c("Data sources", "Preprocessing", "Analysis branches", "Outputs")
lane_x <- c(0.135, 0.365, 0.620, 0.885)
lane_w <- c(0.235, 0.235, 0.270, 0.235)

wrap_text <- function(text, width = 30) {
  if (!nzchar(text)) return(character())
  strwrap(text, width = width)
}

node_by_id <- function(id) {
  nodes[match(id, nodes$id), ]
}

draw_node <- function(node) {
  x <- node$x
  y <- node$y
  w <- node$w
  h <- node$h
  grid.roundrect(
    x = unit(x, "npc"), y = unit(y, "npc"),
    width = unit(w, "npc"), height = unit(h, "npc"),
    r = unit(0.025, "snpc"),
    gp = gpar(fill = node$fill, col = node$border, lwd = 1.2)
  )

  title_lines <- wrap_text(node$title, width = 27)
  summary_lines <- wrap_text(node$summary, width = 33)
  y_top <- y + h / 2 - 0.027

  grid.text(
    title_lines,
    x = unit(x - w / 2 + 0.014, "npc"),
    y = unit(y_top, "npc"),
    just = c("left", "top"),
    gp = gpar(
      fontfamily = "Helvetica", fontface = "bold", fontsize = 8.5,
      col = colors$ink, lineheight = 0.95
    )
  )

  if (length(summary_lines) > 0) {
    title_height <- 0.020 * length(title_lines)
    grid.text(
      summary_lines,
      x = unit(x - w / 2 + 0.014, "npc"),
      y = unit(y_top - title_height - 0.010, "npc"),
      just = c("left", "top"),
      gp = gpar(
        fontfamily = "Helvetica", fontsize = 7.3,
        col = colors$muted, lineheight = 0.95
      )
    )
  }
}

draw_arrow <- function(from, to) {
  src <- node_by_id(from)
  dst <- node_by_id(to)

  x0 <- src$x + src$w / 2
  y0 <- src$y
  x1 <- dst$x - dst$w / 2
  y1 <- dst$y

  if (dst$x <= src$x) {
    x0 <- src$x
    x1 <- dst$x
  }

  grid.segments(
    x0 = unit(x0, "npc"), y0 = unit(y0, "npc"),
    x1 = unit(x1 - 0.010, "npc"), y1 = unit(y1, "npc"),
    arrow = arrow(type = "closed", length = unit(0.08, "inches")),
    gp = gpar(col = colors$line, lwd = 0.8)
  )
}

pdf(output, width = 11, height = 6.8, family = "Helvetica", useDingbats = FALSE)
grid.newpage()

grid.rect(gp = gpar(fill = "white", col = NA))

grid.text(
  "FlyWire Visual System Connectome Analysis Workflow",
  x = unit(0.045, "npc"), y = unit(0.955, "npc"),
  just = c("left", "top"),
  gp = gpar(fontfamily = "Helvetica", fontface = "bold", fontsize = 15, col = colors$ink)
)
grid.text(
  "High-level manuscript workflow from connectome inputs to figures and statistical summaries",
  x = unit(0.045, "npc"), y = unit(0.915, "npc"),
  just = c("left", "top"),
  gp = gpar(fontfamily = "Helvetica", fontsize = 9, col = colors$muted)
)

for (i in seq_along(lane_labels)) {
  grid.roundrect(
    x = unit(lane_x[i], "npc"), y = unit(0.470, "npc"),
    width = unit(lane_w[i], "npc"), height = unit(0.805, "npc"),
    r = unit(0.015, "snpc"),
    gp = gpar(fill = colors$lane, col = NA)
  )
  grid.text(
    lane_labels[i],
    x = unit(lane_x[i] - lane_w[i] / 2 + 0.010, "npc"),
    y = unit(0.865, "npc"),
    just = c("left", "top"),
    gp = gpar(fontfamily = "Helvetica", fontface = "bold", fontsize = 8.8, col = colors$muted)
  )
}

invisible(mapply(draw_arrow, edges$from, edges$to))
invisible(lapply(seq_len(nrow(nodes)), function(i) draw_node(as.list(nodes[i, ]))))

grid.text(
  "Chart intentionally omits job-level scheduling, file paths, batching, and parameter details.",
  x = unit(0.955, "npc"), y = unit(0.030, "npc"),
  just = c("right", "bottom"),
  gp = gpar(fontfamily = "Helvetica", fontsize = 7, col = colors$muted)
)

dev.off()

message("Wrote workflow chart: ", output)
