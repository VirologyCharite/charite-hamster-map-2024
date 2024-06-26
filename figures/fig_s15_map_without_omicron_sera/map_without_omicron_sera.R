
rm(list = ls())

library(tidyverse)
library(Racmacs)

map <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

make_map <- function(mapsSubset) {
  mapsSubset <- optimizeMap(mapsSubset, 2, 1000, options = list(ignore_disconnected = TRUE))
  mapsSubset <- realignMap(mapsSubset, map)

  mapsSubsetP <- procrustesMap(mapsSubset, map, sera = F)
  
  mapsSubsetTriangulation <- triangulationBlobs(mapsSubset, stress_lim = 1, grid_spacing = 0.05)
  
  mapsSubsetBs <- bootstrapMap(
    mapsSubset,
    "resample",
    bootstrap_repeats = 500,
    bootstrap_ags = TRUE,
    bootstrap_sr = TRUE,
    reoptimize = TRUE,
    optimizations_per_repeat = 500,
    ag_noise_sd = 0.4,
    titer_noise_sd = 0.62,
    options = list()
  )
  
  mapsSubsetBsBlobs <- bootstrapBlobs(mapsSubsetBs, conf.level = 0.68, smoothing = 2, gridspacing = 0.05)
  
  result <- list(
    `p` = mapsSubsetP,
    `triang` = mapsSubsetTriangulation,
    `bs` = mapsSubsetBsBlobs
  )
}

# BA.5 sera only
mapBA5Sr <- removeSera(map, sera = srNames(map)[srNames(map) %in% c("1.1", "1.2", "1.3", "5.1", "5.2", "5.3", "6.1", "6.2", "6.3", "12SE0030", "12SE0031", "12SE0032")])
mapBA5Srnmeasured <- colSums(titerTable(mapBA5Sr) != "*")
mapBA5Srndetectable <- colSums(titerTable(mapBA5Sr) != "*" & !grepl("<", titerTable(mapBA5Sr)))
mapBA5SrnmeasuredAg <- rowSums(titerTable(mapBA5Sr) != "*")
mapBA5SrndetectableAg <- rowSums(titerTable(mapBA5Sr) != "*" & !grepl("<", titerTable(mapBA5Sr)))
mapBA5Sr <- subsetMap(
  mapBA5Sr, 
  sera = (mapBA5Srnmeasured > 3 & mapBA5Srndetectable > 2),
  antigens = (mapBA5SrnmeasuredAg > 3 & mapBA5SrndetectableAg > 2)
)
BA5result <- make_map(mapBA5Sr)

# BA.1 sera only
mapBA1Sr <- removeSera(map, sera = srNames(map)[srNames(map) %in% c("1.1", "1.2", "1.3", "6.1", "6.2", "6.3", "9.1", "9.2", "9.3", "12SE0030", "12SE0031", "12SE0032")])
mapBA1Srnmeasured <- colSums(titerTable(mapBA1Sr) != "*")
mapBA1Srndetectable <- colSums(titerTable(mapBA1Sr) != "*" & !grepl("<", titerTable(mapBA1Sr)))
mapBA1SrnmeasuredAg <- rowSums(titerTable(mapBA1Sr) != "*")
mapBA1SrndetectableAg <- rowSums(titerTable(mapBA1Sr) != "*" & !grepl("<", titerTable(mapBA1Sr)))
mapBA1Sr <- subsetMap(
  mapBA1Sr, 
  sera = (mapBA1Srnmeasured > 3 & mapBA1Srndetectable > 2),
  antigens = (mapBA1SrnmeasuredAg > 3 & mapBA1SrndetectableAg > 2)
)
BA1result <- make_map(mapBA1Sr)


# BA.2 sera only
mapBA2Sr <- removeSera(map, sera = srNames(map)[srNames(map) %in% c("5.1", "5.2", "5.3", "9.1", "9.2", "9.3", "12SE0030", "12SE0031", "12SE0032")])
mapBA2Srnmeasured <- colSums(titerTable(mapBA2Sr) != "*")
mapBA2Srndetectable <- colSums(titerTable(mapBA2Sr) != "*" & !grepl("<", titerTable(mapBA2Sr)))
mapBA2SrnmeasuredAg <- rowSums(titerTable(mapBA2Sr) != "*")
mapBA2SrndetectableAg <- rowSums(titerTable(mapBA2Sr) != "*" & !grepl("<", titerTable(mapBA2Sr)))
mapBA2Sr <- subsetMap(
  mapBA2Sr, 
  sera = (mapBA2Srnmeasured > 3 & mapBA2Srndetectable > 2),
  antigens = (mapBA2SrnmeasuredAg > 3 & mapBA2SrndetectableAg > 2)
)
BA2result <- make_map(mapBA2Sr)


# XBB.2 sera only
mapXBB2Sr <- removeSera(map, sera = srNames(map)[srNames(map) %in% c("5.1", "5.2", "5.3", "9.1", "9.2", "9.3", "1.1", "1.2", "1.3", "6.1", "6.2", "6.3")])
mapXBB2Srnmeasured <- colSums(titerTable(mapXBB2Sr) != "*")
mapXBB2Srndetectable <- colSums(titerTable(mapXBB2Sr) != "*" & !grepl("<", titerTable(mapXBB2Sr)))
mapXBB2SrnmeasuredAg <- rowSums(titerTable(mapXBB2Sr) != "*")
mapXBB2SrndetectableAg <- rowSums(titerTable(mapXBB2Sr) != "*" & !grepl("<", titerTable(mapXBB2Sr)))
mapXBB2Sr <- subsetMap(
  mapXBB2Sr, 
  sera = (mapXBB2Srnmeasured > 3 & mapXBB2Srndetectable > 2),
  antigens = (mapXBB2SrnmeasuredAg > 3 & mapXBB2SrndetectableAg > 2)
)
XBB2result <- make_map(mapXBB2Sr)


# No omicron sera
mapNoOmiSr <- removeSera(map, sera = srNames(map)[srNames(map) %in% c("1.1", "1.2", "1.3", "6.1", "6.2", "6.3", "5.1", "5.2", "5.3", "9.1", "9.2", "9.3", "12SE0030", "12SE0031", "12SE0032")])
mapNoOmiSrnmeasured <- colSums(titerTable(mapNoOmiSr) != "*")
mapNoOmiSrndetectable <- colSums(titerTable(mapNoOmiSr) != "*" & !grepl("<", titerTable(mapNoOmiSr)))
mapNoOmiSrnmeasuredAg <- rowSums(titerTable(mapNoOmiSr) != "*")
mapNoOmiSrndetectableAg <- rowSums(titerTable(mapNoOmiSr) != "*" & !grepl("<", titerTable(mapNoOmiSr)))
mapNoOmiSr <- subsetMap(
  mapNoOmiSr, 
  sera = (mapNoOmiSrnmeasured > 3 & mapNoOmiSrndetectable > 2),
  antigens = (mapBA2SrnmeasuredAg > 3 & mapNoOmiSrndetectableAg > 2)
)
noOmiSrresult <- make_map(mapNoOmiSr)

# Plot figure

xmin <- -5
xmax <- 7
ymin <- -3
ymax <- 5
png("figures/fig_s15_map_without_omicron_sera/map_without_omicron_sera_240123.png", width = 15, height = 17.6, units = "in", res=300, pointsize = 18)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), ncol=3))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
plot(noOmiSrresult$p, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-5, 4.4, "A   no Omicron sera", cex = 2, pos = 4)
plot(BA1result$p, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-5, 4.4, "D   no BA.2, BA.5, XBB.2 sera", cex = 2, pos = 4)
plot(BA2result$p, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-5, 4.4, "G   no BA.1, BA.5, XBB.2 sera", cex = 2, pos = 4)
plot(BA5result$p, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-5, 4.4, "J   no BA.1, BA.2, XBB.2 sera", cex = 2, pos = 4)
plot(XBB2result$p, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-5, 4.4, "M   no BA.1, BA.2, BA.5 sera", cex = 2, pos = 4)

plot(noOmiSrresult$triang, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-4.5, 4.5, "B", cex = 2)
plot(BA1result$triang, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-4.5, 4.5, "E", cex = 2)
plot(BA2result$triang, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-4.5, 4.5, "H", cex = 2)
plot(BA5result$triang, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-4.5, 4.5, "K", cex = 2)
plot(XBB2result$triang, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-4.5, 4.5, "N", cex = 2)

plot(noOmiSrresult$bs, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-4.5, 4.5, "C", cex = 2)
plot(BA1result$bs, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-4.5, 4.5, "F", cex = 2)
plot(BA2result$bs, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-4.5, 4.5, "I", cex = 2)
plot(BA5result$bs, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-4.5, 4.5, "L", cex = 2)
plot(XBB2result$bs, xlim = c(xmin, xmax), ylim = c(ymin, ymax), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.5)
text(-4.5, 4.5, "O", cex = 2)
dev.off()

