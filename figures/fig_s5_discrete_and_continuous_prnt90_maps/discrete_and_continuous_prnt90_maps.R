
rm(list = ls())

library(Racmacs)

mapDisc <- rotateMap(read.acmap("data/240123-hamster/map-discrete-90.ace"), 0)
mapContFixTopFixBottom <- rotateMap(read.acmap("data/240123-hamster/map-continuous-fixtop-fixbottom-90.ace"), 0)
mapContFixTop <- rotateMap(read.acmap("data/240123-hamster/map-continuous-fixtop-90.ace"), 0)
mapContFixBottom <- rotateMap(read.acmap("data/240123-hamster/map-continuous-fixbottom-90.ace"), 0)
mapCont <- rotateMap(read.acmap("data/240123-hamster/map-continuous-90.ace"), 0)

sera = T
pDisc <- procrustesMap(mapDisc, mapDisc, sera = sera)
pContFixTopFixBottom <- procrustesMap(mapContFixTopFixBottom, mapDisc, sera = sera)
pContFixTop <- procrustesMap(mapContFixTop, mapDisc, sera = sera)
pContFixBottom <- procrustesMap(mapContFixBottom, mapDisc, sera = sera)
pCont <- procrustesMap(mapCont, mapDisc, sera = sera)

xmin <- -4
xmax <- 8
ymin <- -3
ymax <- 5
png("figures/fig_s5_discrete_and_continuous_prnt90_maps/discrete_and_continuous_prnt90_maps_update_240123.png", width = 14, height = 5.2, units = "in", res=300, pointsize = 18)
layout(matrix(c(1, 5, 2, 0, 3, 0, 4, 0), ncol=4))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
plot(pContFixTopFixBottom, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "A", cex = 2)
text(3, -2.5, "Continuous, fixtop, fixbottom", cex = 1)

plot(pContFixTop, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "B", cex = 2)
text(3.5, -2.5, "Continuous, fixtop", cex = 1)

plot(pContFixBottom, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "C", cex = 2)
text(3.5, -2.5, "Continuous, fixbottom", cex = 1)

plot(pCont, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "D", cex = 2)
text(3.5, -2.5, "Continuous", cex = 1)

plot(pDisc, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5, plot_stress = T)
text(-3.5, 4.5, "E", cex = 2)
text(3.5, -2.5, "Discrete", cex = 1)

dev.off()

# Get RMSDs
procrustesData(mapContFixTopFixBottom, mapDisc, sera = T)$total_rmsd
procrustesData(mapContFixTop, mapDisc, sera = T)$total_rmsd
procrustesData(mapContFixBottom, mapDisc, sera = T)$total_rmsd
procrustesData(mapCont, mapDisc, sera = T)$total_rmsd

procrustesData(mapContFixTopFixBottom, mapContFixTop, sera = T)$total_rmsd
procrustesData(mapContFixTopFixBottom, mapContFixBottom, sera = T)$total_rmsd
procrustesData(mapContFixTopFixBottom, mapCont, sera = T)$total_rmsd

procrustesData(mapContFixTop, mapContFixBottom, sera = T)$total_rmsd
procrustesData(mapContFixTop, mapCont, sera = T)$total_rmsd

procrustesData(mapContFixBottom, mapCont, sera = T)$total_rmsd

