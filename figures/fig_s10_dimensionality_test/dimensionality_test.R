
rm(list = ls())

library(tidyverse)
library(Racmacs)

source("Rcivaclib/viewer_plotting.R")

map <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")


# Make the 3D map
map3D <- optimizeMap(map, 3, 500)
map3Dp <- procrustesMap(map3D, map, sera = F)

map3d_viewer <- plot_map(
  map = procrustes_2d_to_3d(map, map3D),
  scaling_size = 0.15,
  sr_scaling_size = 0.1,
  alter_opacity = F,
  sera_opacity = 0.2,
  grid.col = NA,
  rotation = c(-0.9954, 0.0373, 0.2036),
  translation = c(-0.0799, 0.0153, 0.0057),
  zoom = 1.4
)

view(map3d_viewer)


# Run the dimensionality testing
mapDims <- dimensionTestMap(
  map,
  dimensions_to_test = 1:5,
  test_proportion = 0.1,
  minimum_column_basis = "none",
  fixed_column_bases = colBases(map),
  number_of_optimizations = 100, # 500,
  replicates_per_dimension = 500, # 1000,
  options = list()
)

# > mapDims
# dimensions mean_rmse_detectable var_rmse_detectable mean_rmse_nondetectable var_rmse_nondetectable replicates
# 1          1             1.691746          0.18252833                2.182538              0.5714953        500
# 2          2             1.130425          0.03835408                1.677822              0.2435180        500
# 3          3             1.082180          0.03249740                1.452295              0.1666108        500
# 4          4             1.101438          0.03473201                1.420816              0.1654450        500
# 5          5             1.114212          0.03538285                1.409631              0.1643328        500

# Plot the figure
df <- data.frame(dimensions=c(1, 2, 3, 4, 5),
                 rmse=c(mapDims$mean_rmse_detectable))

mapDims %>%
  ggplot(
    aes(
      x=dimensions,
      y=mean_rmse_detectable,
      ymin=mean_rmse_detectable-var_rmse_detectable,
      ymax=mean_rmse_detectable+var_rmse_detectable)
  ) +
  geom_pointrange() +
  theme_bw() +
  coord_cartesian(ylim = c(1, 2)) +
  labs(
    y="Mean RMSE of detectable titers",
    x="Dimensions"
  )

ggsave("figures/fig_s10_dimensionality_test/dimensionality_test_240123.png", width = 4, height = 3)


# Print the results
print(mapDims$mean_rmse_detectable)
# 1.554070 1.110754 1.069089 1.056519 1.059710
