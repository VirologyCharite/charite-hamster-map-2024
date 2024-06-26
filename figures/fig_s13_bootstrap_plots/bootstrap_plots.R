
rm(list = ls())

library(Racmacs)

neut <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

### Noisy bootstrap
# Titer and antigen noise
neutBootTA <- bootstrapMap(
  neut,
  "noisy",
  bootstrap_repeats = 1000,
  bootstrap_ags = TRUE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = 500,
  ag_noise_sd = 0.4,
  titer_noise_sd = 0.33,
  options = list()
)

neutBootTABlobs <- bootstrapBlobs(neutBootTA, conf.level = 0.68, smoothing = 2, gridspacing = 0.05)

# Titer noise
neutBootT <- bootstrapMap(
  neut,
  "noisy",
  bootstrap_repeats = 1000,
  bootstrap_ags = TRUE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = 500,
  ag_noise_sd = 0,
  titer_noise_sd = 0.33,
  options = list()
)

neutBootTBlobs <- bootstrapBlobs(neutBootT, conf.level = 0.68, smoothing = 2, gridspacing = 0.05)

# Antigen noise
neutBootA <- bootstrapMap(
  neut,
  "noisy",
  bootstrap_repeats = 1000,
  bootstrap_ags = TRUE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = 500,
  ag_noise_sd = 0.4,
  titer_noise_sd = 0,
  options = list()
)

neutBootABlobs <- bootstrapBlobs(neutBootA, conf.level = 0.68, smoothing = 2, gridspacing = 0.05)


## Resample bootstrap

# Resampling antigens and sera
neutBootASR <- bootstrapMap(
  neut,
  "resample",
  bootstrap_repeats = 1000,
  bootstrap_ags = TRUE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = 500,
  options = list()
)

neutBootASBlobsR <- bootstrapBlobs(neutBootASR, conf.level = 0.68, smoothing = 2, gridspacing = 0.05)


# Resampling antigens
neutBootAR <- bootstrapMap(
  neut,
  "resample",
  bootstrap_repeats = 1000,
  bootstrap_ags = TRUE,
  bootstrap_sr = FALSE,
  reoptimize = TRUE,
  optimizations_per_repeat = 500,
  options = list()
)

neutBootABlobsR <- bootstrapBlobs(neutBootAR, conf.level = 0.68, smoothing = 2, gridspacing = 0.05)


# Resampling sera
neutBootSR <- bootstrapMap(
  neut,
  "resample",
  bootstrap_repeats = 1000,
  bootstrap_ags = FALSE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = 500,
  options = list()
)

neutBootSBlobsR <- bootstrapBlobs(neutBootSR, conf.level = 0.68, smoothing = 2, gridspacing = 0.05)




# Set plotting limits
plotlims <- Racmacs:::mapPlotLims(neutBootTABlobs, sera = F)

# Plot the figure
png("figures/fig_s13_bootstrap_plots/bootstrap_plots_240123.png", width = 10, height = 5, units = "in", res=300, pointsize = 18)
layout(matrix(c(1, 4, 2, 5, 3, 6), ncol=3))
#par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
par(oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0))
plot(neutBootTABlobs, xlim = plotlims$xlim, ylim = plotlims$ylim + c(-1, 1), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(plotlims$xlim[1]+0.5, plotlims$ylim[2]+0.5, "A", cex = 2)
plot(neutBootTBlobs, xlim = plotlims$xlim, ylim = plotlims$ylim + c(-1, 1), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(plotlims$xlim[1]+0.5, plotlims$ylim[2]+0.5, "B", cex = 2)
plot(neutBootABlobs, xlim = plotlims$xlim, ylim = plotlims$ylim + c(-1, 1), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(plotlims$xlim[1]+0.5, plotlims$ylim[2]+0.5, "C", cex = 2)


plot(neutBootASBlobsR, xlim = plotlims$xlim, ylim = plotlims$ylim + c(-1, 1), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(plotlims$xlim[1]+0.5, plotlims$ylim[2]+0.5, "D", cex = 2)
plot(neutBootABlobsR, xlim = plotlims$xlim, ylim = plotlims$ylim + c(-1, 1), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(plotlims$xlim[1]+0.5, plotlims$ylim[2]+0.5, "E", cex = 2)
plot(neutBootSBlobsR, xlim = plotlims$xlim, ylim = plotlims$ylim + c(-1, 1), fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(plotlims$xlim[1]+0.5, plotlims$ylim[2]+0.5, "F", cex = 2)
dev.off()

