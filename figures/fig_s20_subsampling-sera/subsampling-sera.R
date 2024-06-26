
rm(list = ls())

library(tidyverse)
library(Racmacs)
library(patchwork)

# Read in map
mapOrig <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

sr_groups_ind <- c("BA.2", "BA.2", "BA.2",
                   "XBB.2", "XBB.2", "XBB.2",
                   "Alpha", "Alpha",
                   "Beta", "Beta", "Beta",
                   "Delta", "Delta", "Delta",
                   "BA.1", "BA.1", "BA.1",
                   "BA.2-12", "BA.2-12", "BA.2-12",
                   "B.1+E484K", "B.1+E484K", "B.1+E484K",
                   "D614G", "D614G", "D614G",
                   "BA.5", "BA.5", "BA.5")

srGroups(mapOrig) <- sr_groups_ind

map <- subsetMap(mapOrig, antigens = c("Alpha", "BA.1", "BA.2-12",
                                       "BA.2", "BA.4", "BA.5", "Beta",
                                       "Delta", "Mu", "D614G", "B.1+E484K"),
                 sera = c("1.1", "1.2", "1.3", "2.1", "2.2", "3.1", "3.2", "3.3",
                          "4.1", "4.2", "4.3", "5.1", "5.2", "5.3", "6.1", "6.2", "6.3",
                          "7.1", "7.2", "7.3", "8.1", "8.2", "8.3"))

map <- optimizeMap(map, 2, 500)

map <- realignMap(map, mapOrig)

# collapse serum groups
srGroups(map) <- forcats::fct_drop(
  forcats::fct_collapse(
    srGroups(map),
    "BA.2" = c("BA.2", "BA.2-12")
  )
)

# Define some functions
plot_remove_srs <- function(map, srs) {
  map %>% removeSera(srs) %>%
    optimizeMap(2, 100, "none", options = list(ignore_disconnected = TRUE)) %>%
    realignMap(map) %>%
    procrustesMap(map, sera = FALSE) -> map_subset
  
  map_subset
  
}


plot_remove_srgroups <- function(map, sr_groups) {
  
  map %>% removeSera(
    !(srGroups(map) %in% sr_groups)
  ) %>%
    optimizeMap(2, 100, "none", options = list(ignore_disconnected = TRUE)) %>%
    realignMap(map) %>%
    procrustesMap(map, sera = FALSE) -> map_subset
  
    map_subset
  
}


# function to randomise a map
randomise_map <- function(map, mapfilename) {
  
  titers <- titerTable(map)
  titers[] <- sample(titers, size = length(titers), replace = FALSE)
  
  titerTable(map) <- titers
  mapNew <- suppressMessages(optimizeMap(map, 2, 100, "none", options = list(ignore_disconnected = TRUE)))
  
  mapNew <- realignMap(mapNew, map)
  
  save.acmap(mapNew, mapfilename)
  
}


## Randomly sample any sera
# Make the subsampled maps
for (srCount in 1:22) {
  for (repeat_ in 1:10) {
    selectedSrs <- sample(srNames(map), srCount)
    
    map_subset <- plot_remove_srs(map, selectedSrs)
    save.acmap(map_subset, paste("figures/fig_s20_subsampling-sera/", srCount, "-", repeat_, ".ace", sep = ""))
  }
}


# Plot the subsampled maps
png("figures/fig_s20_subsampling-sera/subsample_sera.png", width = 34, height = 70, units = "in", res=300, pointsize = 18)
layout(matrix(c(1:210), ncol=10, nrow=21, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
for (srCount in 1:21) {
  print(srCount)
  for (repeat_ in 1:10) {
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/", srCount, "-", repeat_, ".ace", sep = ""))
    map_subsetP <- procrustesMap(map_subset, map, sera = F)
    
    plot(map_subsetP, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.4)
    text(-3.95, 4.5, paste("N sera: ", numSera(map_subset), sep = ""), pos = 4, cex = 1)
  }
}

dev.off()

# Make randomised maps
for (srCount in 1:22) {
  print(srCount)
  for (repeat_ in 1:10) {
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/", srCount, "-", repeat_, ".ace", sep = ""))
    randomise_map(map_subset, paste("figures/fig_s20_subsampling-sera/randomised-", srCount, "-", repeat_, ".ace", sep = ""))
  }
}

# plot rmsd vs number of sera
sr_rmsd <- tibble(
  num_sr = numeric(0),
  rmsd = numeric(0),
  color = character(0),
  type = character(0)
)

for (srCount in 1:22) {
  print(srCount)
  for (repeat_ in 1:10) {
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/", srCount, "-", repeat_, ".ace", sep = ""))
    total_rmsd <- procrustesData(map_subset, map, sera = F)$total_rmsd
    s <- srGroups(map_subset)
    if (("BA.1" %in% s) & ("BA.2" %in% s)) {
      c = "w_omi"
    } else {
      c = "wo_omi"
    }
    
    sr_rmsd <- bind_rows(
      sr_rmsd,
      tibble(
        num_sr = numSera(map_subset),
        rmsd = total_rmsd,
        color = c,
        type = "normal"
        )
      )
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/randomised-", srCount, "-", repeat_, ".ace", sep = ""))
    total_rmsd <- procrustesData(map_subset, map, sera = F)$total_rmsd
    s <- srGroups(map_subset)
    if (("BA.1" %in% s) & ("BA.2" %in% s)) {
      c = "w_omi"
    } else {
      c = "wo_omi"
    }
    
    sr_rmsd <- bind_rows(
      sr_rmsd,
      tibble(
        num_sr = numSera(map_subset),
        rmsd = total_rmsd,
        color = c,
        type = "random"
      )
    )
    
  }
}


# Randomly sample sera within serum group

# Make subsampled maps
srBySrGroup <- list(c("1.1", "1.2", "1.3", "6.1", "6.2", "6.3"), c("2.1", "2.2"), c("3.1", "3.2", "3.3"),
                    c("4.1", "4.2", "4.3"), c("5.1", "5.2", "5.3"),
                    c("7.1", "7.2", "7.3"), c("8.1", "8.2", "8.3"))

for (srCount in 1:2) {
  for (repeat_ in 1:10) {
    
    selectedSrs <- c()
    for (s in srBySrGroup) {
      selectedSr <- sample(s, srCount, replace = FALSE)
      selectedSrs <- c(selectedSrs, c(selectedSr))
    
    positiveSelectedSera <- setdiff(srNames(map), selectedSrs)
    map_subset <- plot_remove_srs(map, positiveSelectedSera)
    save.acmap(map_subset, paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-", srCount, "-", repeat_, ".ace", sep = ""))
    }
  }
}


# Plot the subsampled maps
png("figures/fig_s20_subsampling-sera/subsample_sera_by_sr_group.png", width = 34, height = 7, units = "in", res=300, pointsize = 18)
layout(matrix(c(1:20), ncol=10, nrow=2, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
for (srCount in 1:2) {
  for (repeat_ in 1:10) {
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-", srCount, "-", repeat_, ".ace", sep = ""))
    map_subsetP <- procrustesMap(map_subset, map, sera = F)
    
    plot(map_subsetP, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.4)
    text(-3.95, 4.5, paste("N sera: ", numSera(map_subset), sep = ""), pos = 4, cex = 1)
  }
}

dev.off()




# Randomly sample serum groups
for (srGCount in 2:6) {
  for (repeat_ in 1:10) {
    selectedSrs <- sample(levels(srGroups(map)), srGCount, replace = FALSE)
    
    map_subset <- plot_remove_srgroups(map, selectedSrs)
    save.acmap(map_subset, paste("figures/fig_s20_subsampling-sera/subsample_all_sr_groups-", srGCount, "-", repeat_, ".ace"))
    randomise_map(map_subset, paste("figures/fig_s20_subsampling-sera/randomised-subsample_all_sr_groups-", srGCount, "-", repeat_, ".ace"))
  }
}


png("figures/fig_s20_subsampling-sera/subsample_serum_groups.png", width = 34, height = 18, units = "in", res=300, pointsize = 18)
layout(matrix(c(1:50), ncol=10, nrow=5, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0), cex.main = 0.5)
for (srGCount in 2:6) {
  for (repeat_ in 1:10) {
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsample_all_sr_groups-", srGCount, "-", repeat_, ".ace"))
    
    plot(map_subset, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.4)
    text(-3.95, 4.5, paste(unique(srGroups(map)), collapse="-"), pos = 4, cex = 1)
  }
}

dev.off()


# Subsample serum groups and subsample to two sera per serum group
# Make some maps
srBySrGroupBA2 <- list(c("1.1", "1.2", "1.3", "6.1", "6.2", "6.3"), c("2.1", "2.2"), c("3.1", "3.2", "3.3"),
                    c("4.1", "4.2", "4.3"), c("5.1", "5.2", "5.3"),
                    c("7.1", "7.2", "7.3"), c("8.1", "8.2", "8.3"))

for (srGCount in 3:7) {
  for (repeat_ in 1:10) {
    
    # Subsample the serum groups
    selectedSrGroups <- sample(srBySrGroupBA2, srGCount, replace = FALSE)
    
    # Subsample the sera to two sr per sr group
    selectedSrs <- c()
    for (s in selectedSrGroups) {
      selectedSr <- sample(s, 2, replace = FALSE)
      selectedSrs <- c(selectedSrs, c(selectedSr))
      
      positiveSelectedSera <- setdiff(srNames(map), selectedSrs)
      
    }
    print(paste(srGCount, repeat_))
    map_subset <- plot_remove_srs(map, positiveSelectedSera)
    save.acmap(map_subset, paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-2sr-", srGCount, "-", repeat_, ".ace", sep = ""))
  }
}

# plot
png("figures/fig_s20_subsampling-sera/subsample_srgrous_2sr.png", width = 40, height = 20, units = "in", res=300, pointsize = 18)
layout(matrix(c(1:50), ncol=10, nrow=5, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
for (srCount in 3:7) {
  for (repeat_ in 1:10) {
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-2sr-", srCount, "-", repeat_, ".ace", sep = ""))
    map_subsetP <- procrustesMap(map_subset, map, sera = F)
    
    plot(map_subsetP, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.4)
    text(-3.95, 4.5, paste("N sera: ", numSera(map_subset), "; N srGroups: ", numSeraGroups(map_subset), sep = ""), pos = 4, cex = 1)
  }
}

dev.off()



# Subsample serum groups and subsample to one serum per serum group

# Make some maps
srBySrGroupBA2 <- list(c("1.1", "1.2", "1.3", "6.1", "6.2", "6.3"), c("2.1", "2.2"), c("3.1", "3.2", "3.3"),
                       c("4.1", "4.2", "4.3"), c("5.1", "5.2", "5.3"),
                       c("7.1", "7.2", "7.3"), c("8.1", "8.2", "8.3"))

for (srGCount in 3:7) {
  for (repeat_ in 1:10) {
    
    # Subsample the serum groups
    selectedSrGroups <- sample(srBySrGroupBA2, srGCount, replace = FALSE)
    
    # Subsample the sera to two sr per sr group
    selectedSrs <- c()
    for (s in selectedSrGroups) {
      selectedSr <- sample(s, 1, replace = FALSE)
      selectedSrs <- c(selectedSrs, c(selectedSr))
      
      positiveSelectedSera <- setdiff(srNames(map), selectedSrs)
      
    }
    map_subset <- plot_remove_srs(map, positiveSelectedSera)
    save.acmap(map_subset, paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-1sr-", srGCount, "-", repeat_, ".ace", sep = ""))
  }
}

png("figures/fig_s20_subsampling-sera/subsample_srgrous_1sr.png", width = 40, height = 20, units = "in", res=300, pointsize = 18)
layout(matrix(c(1:50), ncol=10, nrow=5, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
for (srCount in 3:7) {
  for (repeat_ in 1:10) {
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-1sr-", srCount, "-", repeat_, ".ace", sep = ""))
    map_subsetP <- procrustesMap(map_subset, map, sera = F)
    
    plot(map_subsetP, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.4)
    text(-3.95, 4.5, paste("N sera: ", numSera(map_subset), "; N srGroups: ", numSeraGroups(map_subset), sep = ""), pos = 4, cex = 1)
  }
}

dev.off()


# Exhaustively sample serum groups
for (count in 1:5) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    srGToKeep <- c(x[, index], c("BA.1", "BA.2"))
    map_subset <- plot_remove_srgroups(map, srGToKeep)
    save.acmap(map_subset, paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-", count, "-", index, ".ace", sep = ""))
  }
}

# Randomise maps
for (count in 1:5) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    map_subset_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-", count, "-", index, ".ace", sep = ""))
    randomise_map(map_subset_subset, paste("figures/fig_s20_subsampling-sera/randomised-subsampled-srgroups-keep-", count, "-", index, ".ace", sep = ""))
  }
}


png("figures/fig_s20_subsampling-sera/subsample_nonOmi_srgrous.png", width = 29, height = 36, units = "in", res=300, pointsize = 18)
layout(matrix(c(c(1:5), c(0, 0, 0, 0, 0), c(6:15), c(0, 0, 0, 0, 0), c(16:25), c(0, 0, 0, 0, 0), c(26:30)), ncol=5, nrow=9, byrow = T), 
       heights=c(3, 0.5, 3, 3, 0.5, 3, 3, 0.5, 3))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
for (count in 1:4) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-", count, "-", index, ".ace", sep = ""))
    map_subsetP <- procrustesMap(map_subset, map, sera = F)
    total_rmsd <- procrustesData(map_subset, map, sera = F)$total_rmsd
    
    #plot(neutBootSBlobs, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.7)
    plot(map_subsetP, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.7)
    text(-3.95, 4.5, paste("Present: ", paste(unique(srGroups(map_subset))[!(unique(srGroups(map_subset)) %in% c("BA.1", "BA.2"))], collapse=", ")), pos = 4, cex = 1.7)
    text(3.5, -2.5, round(total_rmsd, 2), cex = 1.7)
  }
}

dev.off()

# Subsample exhaustively subsampled serum groups down to one serum per group
srGInfo <- list(
  "BA.2" = c("1.1", "1.2", "1.3", "6.1", "6.2", "6.3"),
  "Alpha" = c("2.1", "2.2"),
  "Beta" = c("3.1", "3.2", "3.3"),
  "Delta" = c("4.1", "4.2", "4.3"),
  "BA.1" = c("5.1", "5.2", "5.3"),
  "B.1+E484K" = c("7.1", "7.2", "7.3"),
  "D614G" = c("8.1", "8.2", "8.3")
)


for (count in 1:5) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    map_subset = read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-", count, "-", index, ".ace", sep = ""))
    # Gather sera information
    ss <- list()
    for (srG in unique(srGroups(map_subset))) {
      ss <- append(ss, srGInfo[srG])
    }
    for (repeat_ in 1:3) {
      selectedSera <- c()
      for (srG in ss) {
        selectedSera <- c(selectedSera, sample(srG, 1))
      }
      
      srToRemove <- srNames(map)[!(srNames(map) %in% selectedSera)]
      map_subset_subset <- plot_remove_srs(map, srToRemove)
      save.acmap(map_subset_subset, paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-1sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
    }
  }
}

# Randomise those maps
for (count in 1:5) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    for (repeat_ in 1:3) {
      map_subset_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-1sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
      randomise_map(map_subset_subset, paste("figures/fig_s20_subsampling-sera/randomised-subsampled-srgroups-keep-1sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
    }
  }
}



# Subsample exhaustively subsampled serum groups down to two sera per group

srGInfo <- list(
  "BA.2" = c("1.1", "1.2", "1.3", "6.1", "6.2", "6.3"),
  "Alpha" = c("2.1", "2.2"),
  "Beta" = c("3.1", "3.2", "3.3"),
  "Delta" = c("4.1", "4.2", "4.3"),
  "BA.1" = c("5.1", "5.2", "5.3"),
  "B.1+E484K" = c("7.1", "7.2", "7.3"),
  "D614G" = c("8.1", "8.2", "8.3")
)


for (count in 1:5) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    map_subset = read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-", count, "-", index, ".ace", sep = ""))
    # Gather sera information
    ss <- list()
    for (srG in unique(srGroups(map_subset))) {
      ss <- append(ss, srGInfo[srG])
    }
    for (repeat_ in 1:3) {
      selectedSera <- c()
      for (srG in ss) {
        selectedSera <- c(selectedSera, sample(srG, 2))
      }
      
      srToRemove <- srNames(map)[!(srNames(map) %in% selectedSera)]
      map_subset_subset <- plot_remove_srs(map, srToRemove)
      save.acmap(map_subset_subset, paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-2sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
    }
  }
}

# randomise those maps
for (count in 1:5) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    for (repeat_ in 1:3) {
      map_subset_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-2sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
      randomise_map(map_subset_subset, paste("figures/fig_s20_subsampling-sera/randomised-subsampled-srgroups-keep-2sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
    }
  }
}

# Plot the exhaustively sampled serum groups subsamples
png("figures/fig_s20_subsampling-sera/subsample_nonOmi_srgrous_subsampled.png", width = 25, height = 110, units = "in", res=300, pointsize = 18)
layout(matrix(1:210, ncol=7, nrow=30, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
for (count in 1:4) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-", count, "-", index, ".ace", sep = ""))
    map_subsetP <- procrustesMap(map_subset, map, sera = F)
    total_rmsd <- procrustesData(map_subset, map, sera = F)$total_rmsd
    
    #plot(neutBootSBlobs, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.7)
    plot(map_subsetP, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.7)
    text(-3.95, 4.5, paste("Present: ", paste(unique(srGroups(map_subset))[!(unique(srGroups(map_subset)) %in% c("BA.1", "BA.2"))], collapse=", ")), pos = 4, cex = 1)
    text(3.5, -2.5, round(total_rmsd, 2), cex = 1.7)
    
    # Plot the subsampled to two sera
    for (repeat_ in 1:3) {
      map_subset_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-2sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
      map_subsetP <- procrustesMap(map_subset_subset, map, sera = F)
      total_rmsd <- procrustesData(map_subset_subset, map, sera = F)$total_rmsd
      
      plot(map_subsetP, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.7)
      text(-3.95, 4.5, paste("2"), pos = 4, cex = 1.7)
      text(3.5, -2.5, round(total_rmsd, 2), cex = 1.7)
    }
    
    # Plot the subsampled to one serum
    for (repeat_ in 1:3) {
      map_subset_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-1sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
      map_subsetP <- procrustesMap(map_subset_subset, map, sera = F)
      total_rmsd <- procrustesData(map_subset_subset, map, sera = F)$total_rmsd
      
      plot(map_subsetP, xlim = c(-4, 4), ylim = c(-3, 5), cex = 0.7)
      text(-3.95, 4.5, paste("1"), pos = 4, cex = 1.7)
      text(3.5, -2.5, round(total_rmsd, 2), cex = 1.7)
    }
  }
}

dev.off()


## BIG RMSD plot

# Gather the data
sr_rmsd <- tibble(
  num_sr = numeric(0),
  rmsd = numeric(0),
  omi = character(0),
  random = character(0),
  type = character(0)
)

# Gather the data from the subsampled individual sera
for (srCount in 1:22) {
  print(srCount)
  for (repeat_ in 1:10) {
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/", srCount, "-", repeat_, ".ace", sep = ""))
    total_rmsd <- procrustesData(map_subset, map, sera = F)$total_rmsd
    s <- srGroups(map_subset)
    if (("BA.1" %in% s) & ("BA.2" %in% s)) {
      c = "w_omi"
    } else {
      c = "wo_omi"
    }
    
    sr_rmsd <- bind_rows(
      sr_rmsd,
      tibble(
        num_sr = numSera(map_subset),
        rmsd = total_rmsd,
        omi = c,
        random = "normal",
        type = "subsample-sera"
      )
    )
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/randomised-", srCount, "-", repeat_, ".ace", sep = ""))
    total_rmsd <- procrustesData(map_subset, map, sera = F)$total_rmsd
    s <- srGroups(map_subset)
    if (("BA.1" %in% s) & ("BA.2" %in% s)) {
      c = "w_omi"
    } else {
      c = "wo_omi"
    }
    
    sr_rmsd <- bind_rows(
      sr_rmsd,
      tibble(
        num_sr = numSera(map_subset),
        rmsd = total_rmsd,
        omi = c,
        random = "random",
        type = "subsample-sera"
      )
    )
  }
}

# Gather the data from the subsampled serum groups, without BA.1 and BA.2 sera
for (count in 1:4) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    
    m <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-", count, "-", index, ".ace", sep = ""))
    total_rmsd <- procrustesData(m, map, sera = F)$total_rmsd
    
    sr_rmsd <- bind_rows(
      sr_rmsd,
      tibble(
        num_sr = numSera(m),
        rmsd = total_rmsd,
        omi = "w_omi",
        random = "normal",
        type = "subsample-sr_groups"
      )
    )
    
    m <- read.acmap(paste("figures/fig_s20_subsampling-sera/randomised-subsampled-srgroups-keep-", count, "-", index, ".ace", sep = ""))
    total_rmsd <- procrustesData(m, map, sera = F)$total_rmsd
    
    sr_rmsd <- bind_rows(
      sr_rmsd,
      tibble(
        num_sr = numSera(m),
        rmsd = total_rmsd,
        omi = "w_omi",
        random = "random",
        type = "subsample-sr_groups"
      )
    )
    
  }
}


# Gather the data from the subsampled serum groups, without BA.1 and BA.2 sera, subsampling to two sera
for (count in 1:5) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    for (repeat_ in 1:3) {
      map_subset_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-2sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
      
      total_rmsd <- procrustesData(map_subset_subset, map, sera = F)$total_rmsd
      
      sr_rmsd <- bind_rows(
        sr_rmsd,
        tibble(
          num_sr = numSera(map_subset_subset),
          rmsd = total_rmsd,
          omi = "w_omi",
          random = "normal",
          type = "subsample-sr_groups-2sr"
        )
      )
      
      map_subset_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/randomised-subsampled-srgroups-keep-2sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
      total_rmsd <- procrustesData(map_subset_subset, map, sera = F)$total_rmsd
      
      sr_rmsd <- bind_rows(
        sr_rmsd,
        tibble(
          num_sr = numSera(map_subset_subset),
          rmsd = total_rmsd,
          omi = "w_omi",
          random = "random",
          type = "subsample-sr_groups-2sr"
        )
      )
    }
  }
}


# Gather the data from the subsampled serum groups, without BA.1 and BA.2 sera, subsampling to one serum
for (count in 1:5) {
  x <- combn(levels(srGroups(map))[!(levels(srGroups(map)) %in% c("BA.1", "BA.2"))], count)
  for (index in 1:ncol(x)) {
    for (repeat_ in 1:3) {
      map_subset_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsampled-srgroups-keep-1sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
      
      total_rmsd <- procrustesData(map_subset_subset, map, sera = F)$total_rmsd
      
      sr_rmsd <- bind_rows(
        sr_rmsd,
        tibble(
          num_sr = numSera(map_subset_subset),
          rmsd = total_rmsd,
          omi = "w_omi",
          random = "normal",
          type = "subsample-sr_groups-1sr"
        )
      )
      
      map_subset_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/randomised-subsampled-srgroups-keep-1sr-", count, "-", index, "-", repeat_, ".ace", sep = ""))
      total_rmsd <- procrustesData(map_subset_subset, map, sera = F)$total_rmsd
      
      sr_rmsd <- bind_rows(
        sr_rmsd,
        tibble(
          num_sr = numSera(map_subset_subset),
          rmsd = total_rmsd,
          omi = "w_omi",
          random = "random",
          type = "subsample-sr_groups-1sr"
        )
      )
    }
  }
}

# Gather data for randomly subsampled serum groups
for (srGCount in 2:6) {
  for (repeat_ in 1:10) {
    selectedSrs <- sample(levels(srGroups(map)), srGCount, replace = FALSE)
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/subsample_all_sr_groups-", srGCount, "-", repeat_, ".ace"))
    total_rmsd <- procrustesData(map_subset, map, sera = F)$total_rmsd
    s <- srGroups(map_subset)
    if (("BA.1" %in% s) & ("BA.2" %in% s)) {
      c = "w_omi"
    } else {
      c = "wo_omi"
    }
    
    sr_rmsd <- bind_rows(
      sr_rmsd,
      tibble(
        num_sr = numSera(map_subset),
        rmsd = total_rmsd,
        omi = c,
        random = "normal",
        type = "subsample-sr_groups_all"
      )
    )
    
    map_subset <- read.acmap(paste("figures/fig_s20_subsampling-sera/randomised-subsample_all_sr_groups-", srGCount, "-", repeat_, ".ace"))
    total_rmsd <- procrustesData(map_subset, map, sera = F)$total_rmsd
    
    sr_rmsd <- bind_rows(
      sr_rmsd,
      tibble(
        num_sr = numSera(map_subset),
        rmsd = total_rmsd,
        omi = "w_omi",
        random = "random",
        type = "subsample-sr_groups_all"
      )
    )
  }
}

sr_rmsd %>%
  mutate(
    type = factor(type, levels = c("subsample-sera", "subsample-sr_groups_all", "subsample-sr_groups", "subsample-sr_groups-2sr", "subsample-sr_groups-1sr"))
  ) -> sr_rmsd

x <- c(
  "subsample-sera" = "A) Subsample sera",
  "subsample-sr_groups_all" = "B) Subsample serum groups",
  "subsample-sr_groups" = "C) Subsample pre-Omicron serum groups, all sera",
  "subsample-sr_groups-2sr" = "D) Subsample pre-Omicron serum groups, two sera",
  "subsample-sr_groups-1sr" = "E) Subsample pre-Omicron serum groups, one serum"
)

sr_rmsd %>%
  subset(
    random == "normal"
  ) %>%
  ggplot(
    aes(
      x = factor(num_sr),
      y = rmsd
    )
  ) + 
  geom_boxplot() +
  geom_point(
    aes(
      color = omi
    ),
    #color = "deepskyblue",
    alpha = 0.7
  ) +
  geom_point(
    data = subset(sr_rmsd, random == "random"),
    aes(
      x = factor(num_sr),
      y = rmsd
    ),
    shape = 4,
    alpha = 0.3
  ) +
  theme_bw() +
  scale_color_discrete(labels=c("With >1 Omicron sera", "Without Omicron sera")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 10)) +
  facet_wrap(~type, ncol = 5, labeller = as_labeller(x)) +
  labs(
    x = "Number of sera",
    y = "Root mean squared deviation of variant positions",
    color=""
  ) 
ggsave("figures/fig_s20_subsampling-sera/plot_240123.png", width = 20, height = 4)



