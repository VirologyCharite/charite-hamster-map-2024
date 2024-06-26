
# Setup workspace
rm(list = ls())
library(tidyverse)
library(patchwork)
library(Racmacs)

source("Rcivaclib/common.R")

# Load map
map <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

# Set serum group order
srGroups(map) <- factor(srGroups(map), sr_group_order)

# Set plot limits
maplims <- Racmacs:::mapPlotLims(map, sera = FALSE)

plot_remove_srgroups <- function(map, sr_groups, plottitle) {
  
  map %>% removeSera(
    srGroups(map) %in% sr_groups
  ) -> map_subset
  map_subset <- optimizeMap(map_subset, 2, 500, "none", options = list(ignore_disconnected = TRUE))
  map_subset <- realignMap(map_subset, map) 
  map_subset <- procrustesMap(map_subset, map, sera = FALSE) 
  
  ggplot(
    map_subset,
    xlim = maplims$xlim,
    ylim = c(maplims$ylim[1], maplims$ylim[2]),
    plot_stress = TRUE
  ) + 
    labs(
      title = plottitle
    ) + 
    theme(legend.position = "none")
  
}


plot_remove_variants <- function(map, variant, plottitle) {
  
  map %>% removeAntigens(
    agNames(map) == variant
  ) -> map_subset
  map_subset <- optimizeMap(map_subset, 2, 500, "none")
  map_subset <- realignMap(map_subset, map) 
  map_subset <- procrustesMap(map_subset, map, sera = FALSE) 
  
  ggplot(
    map_subset,
    xlim = maplims$xlim,
    ylim = maplims$ylim,
    plot_stress = TRUE,
    cex = 0.5
  ) + 
    labs(
      title = plottitle
    ) + theme(legend.position = "none")
  
}

# Make maps subsampling sera
maps_sera <- lapply(
  levels(srGroups(map)), function(sr_group) {
    
    # Remove sera from serum group
    message(sr_group)
    plot_remove_srgroups(map, sr_group, paste("Without", sr_group, ""))
    
  }
)

# Remove BA.2 sera
map %>% removeSera(
  c("1.1", "1.2", "1.3", "6.1", "6.2", "6.3")
) %>%
  optimizeMap(2, 500, "none") %>%
  realignMap(map) %>%
  procrustesMap(map, sera = FALSE) -> map_subset

gpA_sera <- ggplot(
  map_subset,
  xlim = maplims$xlim,
  ylim = c(maplims$ylim[1], maplims$ylim[2]),
  plot_stress = TRUE,
  cex = 0.1
) + 
  labs(
    title = "w/o all BA.2 sera"
  )


# Make maps subsampling variants
# Set plot limits
maplims <- Racmacs:::mapPlotLims(map, sera = FALSE)

maps_antigen <- lapply(
  c("D614G", "Alpha", "Delta", "Beta", "Mu", "B.1+E484K", "BA.2-12", "BA.2", "BA.1", "BA.4", "BA.5",
    "BF.7", "BQ.1.18", "BN.1.3.1", "XBB.2", "EG.5.1", "JN.1"), function(variant) {
      
      # Remove sera from serum group
      message(variant)
      plot_remove_variants(map, variant, paste("Without", variant, "variant"))
      
    }
)

# Make plots without both BA.2 variants
# Remove BA.2 sera
map %>% removeAntigens(
  c("BA.2", "BA.2-12")
) %>%
  optimizeMap(2, 500, "none") %>%
  realignMap(map) %>%
  procrustesMap(map, sera = FALSE) -> map_subset

gpA_antigen <- ggplot(
  map_subset,
  xlim = maplims$xlim,
  ylim = c(maplims$ylim[1], maplims$ylim[2]),
  plot_stress = TRUE
) + 
  labs(
    title = "w/o all BA.2 variants"
  )


# Do the plotting
design <- "
  ABCD
  EFGJ
  HIK#
  LMNO
  PQRS
  TUVW
  XYZ1
  23##
"

gp_sera <- wrap_plots(c(maps_sera, list(gpA_sera)), 
                 ncol = 4, guides = "collect") + theme(legend.position = "none")

gp_sera <- gp_sera + guides(shape = "none")

gp_sera

ggsave("figures/fig_s14_subsampling_antigens_and_sera/subsampling_antigens_and_sera_sera_240123.png", gp_sera, width = 15, height = 8, units = "in")

gp_antigen <- wrap_plots(c(maps_antigen, list(gpA_antigen)), 
                      ncol = 4, guides = "collect") + theme(legend.position = "none")
gp_antigen <- gp_antigen + guides(shape = "none")

gp_antigen

ggsave("figures/fig_s14_subsampling_antigens_and_sera/subsampling_antigens_and_sera_antigen_240123.png", gp_antigen, width = 15, height = 13, units = "in")


