
rm(list = ls())

library(Racmacs)

map <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

m_ <- rotateMap(map, -20)
m_ <- translateMap(m_, c(0, -0.5))

ptDrawingOrder(m_) <- c(
  46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26,
  25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 2, 4, 3, 5, 1)

cex_ = 0.6

png("figures/fig_2_main_map/main_map_update_240123.png", width = 6.5, height = 4.5, units = "in", res=300,
    pointsize = 18)
par(mar = c(0.3, 0, 0.3, 0))
plot(m_, xlim = c(-4, 8), ylim = c(-3, 5), cex = cex_, fill.alpha = 1, plot_labels = F,
     outline.alpha = 1, grid.col = "#cfcfcf", plot_stress = F,
     grid.margin.col = "black", show_error_lines = F)
dev.off()

# Get euclidean distances

euclidean_dist <- function(x, y){
  sqrt(sum((x - y)^2))
}

coords <- agCoords(map)
# calculate the distances
euclidean_dist(coords["D614G", ], coords["BA.1", ])
euclidean_dist(coords["D614G", ], coords["BA.2", ])
euclidean_dist(coords["D614G", ], coords["BA.4", ])
euclidean_dist(coords["D614G", ], coords["BA.5", ])
euclidean_dist(coords["BA.1", ], coords["BA.4", ])
euclidean_dist(coords["BA.1", ], coords["BA.5", ])

# Get euclidean distances with confidence intervals:

calculate_euclidean_distances <- function(map, ag1, ag2) {
  
  result <- c()
  
  for (opt in 1:100) {
    coords <- agCoords(map, optimization_number = opt)
    dist <- euclidean_dist(coords[ag1, ], coords[ag2, ])
    result <- append(result, dist)
  }
  
  result
  
}

CI(calculate_euclidean_distances(map, "D614G", "BA.1"), ci = 0.95)
CI(calculate_euclidean_distances(map, "D614G", "BA.5"), ci = 0.95)
CI(calculate_euclidean_distances(map, "D614G", "BA.4"), ci = 0.95)
CI(calculate_euclidean_distances(map, "BA.5", "BA.1"), ci = 0.95)
CI(calculate_euclidean_distances(map, "BA.4", "BA.1"), ci = 0.95)


map3D <- optimizeMap(map, 3, 500)
CI(calculate_euclidean_distances(map3D, "D614G", "BA.1"), ci = 0.95)
CI(calculate_euclidean_distances(map3D, "D614G", "BA.5"), ci = 0.95)
CI(calculate_euclidean_distances(map3D, "D614G", "BA.4"), ci = 0.95)
CI(calculate_euclidean_distances(map3D, "BA.5", "BA.1"), ci = 0.95)
CI(calculate_euclidean_distances(map3D, "BA.4", "BA.1"), ci = 0.95)


CI(calculate_euclidean_distances(map, "D614G", "BA.5"), ci = 0.95)

