
rm(list = ls())

library(Racmacs)

# Read in discrete maps
mapd50 <- read.acmap("data/240123-hamster/map-discrete-50.ace")
mapd75 <- read.acmap("data/240123-hamster/map-discrete-75.ace")
mapd90 <- read.acmap("data/240123-hamster/map-discrete-90.ace")
mapd99 <- read.acmap("data/240123-hamster/map-discrete-99.ace")

sera = T
pd50 <- procrustesMap(mapd50, mapd50, sera = sera)
pd75 <- procrustesMap(mapd75, mapd50, sera = sera)
pd90 <- procrustesMap(mapd90, mapd50, sera = sera)
pd99 <- procrustesMap(mapd99, mapd50, sera = sera)

# Read in continuous maps
mapc50 <- read.acmap("data/240123-hamster/map-continuous-fixtop-fixbottom-50.ace")
mapc75 <- read.acmap("data/240123-hamster/map-continuous-fixtop-fixbottom-75.ace")
mapc90 <- read.acmap("data/240123-hamster/map-continuous-fixtop-fixbottom-90.ace")
mapc99 <- read.acmap("data/240123-hamster/map-continuous-fixtop-fixbottom-99.ace")

sera = T
pc50 <- procrustesMap(mapc50, mapc50, sera = sera)
pc75 <- procrustesMap(mapc75, mapc50, sera = sera)
pc90 <- procrustesMap(mapc90, mapc50, sera = sera)
pc99 <- procrustesMap(mapc99, mapc50, sera = sera)



xmin <- -4
xmax <- 8
ymin <- -3
ymax <- 5
png("figures/fig_s3_sensitivity_levels_antigenic_maps/sensitivity_levels_antigenic_maps_discrete_continuous_240123.png",
    width = 6, height = 8, units = "in", res=300, pointsize = 18)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol=2))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
# Discrete maps
plot(pd50, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "A", cex = 1.5)
text(-3, 4.3, "Discrete, NT50", cex = 1, pos = 4)

plot(pd75, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "B", cex = 1.5)
text(-3, 4.3, "Discrete, NT75", cex = 1, pos = 4)

plot(pd90, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "C", cex = 1.5)
text(-3, 4.3, "Discrete, NT90", cex = 1, pos = 4)

plot(pd99, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "D", cex = 1.5)
text(-3, 4.3, "Discrete, NT99", cex = 1, pos = 4)

# Continuous maps
plot(pc50, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "E", cex = 1.5)
text(-3, 4.3, "Continuous, NT50", cex = 1, pos = 4)

plot(pc75, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "F", cex = 1.5)
text(-3, 4.3, "Continuous, NT75", cex = 1, pos = 4)

plot(pc90, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "G", cex = 1.5)
text(-3, 4.3, "Continuous, NT90", cex = 1, pos = 4)

plot(pc99, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "H", cex = 1.5)
text(-3, 4.3, "Continuous, NT99", cex = 1, pos = 4)

dev.off()
