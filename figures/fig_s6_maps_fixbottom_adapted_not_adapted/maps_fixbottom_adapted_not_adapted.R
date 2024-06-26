
rm(list = ls())

library(Racmacs)

map90 <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90.ace")
map90adapted <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

sera = T
p90adapted <- procrustesMap(map90adapted, map90, sera = sera)


xmin <- -4
xmax <- 8
ymin <- -3
ymax <- 5
png("figures/fig_s6_maps_fixbottom_adapted_not_adapted/maps_fixbottom_adapted_not_adapted_update_240123.png", width = 13.5, height = 3.75, units = "in", res=300, pointsize = 18)
layout(matrix(c(1, 2, 3), ncol=3))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
plot(map90, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7, plot_stress = T)
text(-3.5, 4.5, "A", cex = 2)

plot(map90adapted, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7, plot_stress = T)
text(-3.5, 4.5, "B", cex = 2)

plot(p90adapted, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7, plot_stress = T)
text(-3.5, 4.5, "C", cex = 2)

dev.off()

# Get RMSDs
p90adaptedData <- procrustesData(map90adapted, map90, sera = T)
p90adaptedData$total_rmsd
