
rm(list = ls())

library(tidyverse)
library(Racmacs)
library(covutils)
library(patchwork)
library(titertools)

source("Rcivaclib/common.R")

map <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

data <- long_map_info(map)

sr_groups_info <- read.csv("data/metadata/sr_groups.csv", stringsAsFactors = FALSE)

scatter_antigens <- function(ag1, ag2, linecolor) {
  
  xlim <- c(0, 10)

  d <- tibble(
    x = logtiterTable(map)[ag1, ],
    y = logtiterTable(map)[ag2, ],
    "sr_name" = names(logtiterTable(map)[ag1, ]),
    "sr_group" = sr_groups_info$Serum.group[match(names(logtiterTable(map)[ag1, ]), sr_groups_info$Serum)]
  )

  dnl <- tibble(
    x = titerTable(map)[ag1, ],
    y = titerTable(map)[ag2, ],
    "sr_name" = names(logtiterTable(map)[ag1, ]),
    "sr_group" = sr_groups_info$Serum.group[match(names(logtiterTable(map)[ag1, ]), sr_groups_info$Serum)]
  )
  
  titer_differences <- titertools::log2diff(
    titers2 = dnl$x,
    titers1 = dnl$y,
    dilution_stepsize = 0,
    ci_method = "HDI",
    mu_prior_mu = 0,
    mu_prior_sigma = 100,
    sigma_prior_alpha = 2,
    sigma_prior_beta = 0.75
  )

  d %>%
    mutate(
      sr_group = factor(sr_group, levels = sr_group_order)
    ) %>%
    ggplot(
    ) +
    geom_point(
      aes(
        x = x,
        y = y,
        color = sr_group
      ),
      size = 3,
      alpha = 0.7
    ) +
    geom_abline(slope=1, intercept = 0, lty = "dashed", color = "grey") +
    geom_abline(slope=1, intercept = -titer_differences["mean", "estimate"], lty = "solid", color = linecolor) +
    geom_ribbon(
      data = tibble(
        x = extendrange(xlim),
        ymin = extendrange(xlim) - titer_differences["mean", "lower"],
        ymax = extendrange(xlim) - titer_differences["mean", "upper"]
      ),
      aes(
        x = x,
        ymin = ymin,
        ymax = ymax
      ),
      alpha = 0.3,
      fill = linecolor
    ) +
    coord_cartesian(
      ylim = c(0, 10),
      xlim = c(0, 10)
    ) +
    scale_y_continuous(
      breaks = seq(0, 10),
      labels = c("<20", "20", "40", "80", "160", "320", "640", "1280", "2560", "5120", ">5120"),
      minor_breaks = NULL
    ) +
    scale_x_continuous(
      breaks = seq(0, 10),
      labels = c("<20", "20", "40", "80", "160", "320", "640", "1280", "2560", "5120", ">5120"),
      minor_breaks = NULL
    ) +
    scale_color_manual(
      name = "Serum group",
      values = sr_group_cols
    ) +
    theme_bw() +
    labs(
      x = ag1,
      y = ag2
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) -> gp
  
  gp
}


ba2s <- scatter_antigens("BA.2-12", "BA.2", "#8F4F81")

ba45s <- scatter_antigens("BA.4", "BA.5", "#F08DA5")

gps <- ba2s + ba45s + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 20))
ggsave("figures/fig_s7_ba2s_and_ba45_titer_comparison/ba2s_and_ba45_titer_comparison_240123.png", gps, width = 9.5, height = 4)

