
rm(list = ls())

library(tidyverse)
library(Racmacs)

map <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

mapTriangulation <- triangulationBlobs(map, stress_lim = 1, grid_spacing = 0.05)


# Plot figure

xmin <- -4
xmax <- 8
ymin <- -3
ymax <- 4
png("figures/fig_s16_error_and_triangulation_blobs/error_and_triangulation_blobs_240123.png", width = 9, height = 3, units = "in", res=300, pointsize = 18)
layout(matrix(c(1, 2), ncol=2))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
plot(map, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, show_error_lines = T)
text(-3.5, 3.5, "A", cex = 1.5)
plot(mapTriangulation, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3)
text(-3.5, 3.5, "B", cex = 1.5)

dev.off()

