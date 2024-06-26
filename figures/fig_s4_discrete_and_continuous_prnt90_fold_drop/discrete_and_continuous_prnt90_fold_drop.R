
rm(list = ls())

library(tidyverse)
library(Racmacs)
library(patchwork)
library(titertools)

source("Rcivaclib/scales.R")
source("Rcivaclib/calculate_fold_change.R")
source("Rcivaclib/map_longinfo.R")

mapDisc <- read.acmap("data/240123-hamster/map-discrete-90.ace")
mapContFixTopFixBottom <- read.acmap("data/240123-hamster/map-continuous-fixtop-fixbottom-90.ace")
mapContFixTop <- read.acmap("data/240123-hamster/map-continuous-fixtop-90.ace")
mapContFixBottom <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90.ace")
mapCont <- read.acmap("data/240123-hamster/map-continuous-90.ace")

sr_groups <- c("Beta sera", "Delta sera", "D614G sera", "Alpha sera",
               "BA.1 sera", "BA.2 sera", "BA.2-12 sera", "B.1+E484K sera",
               "BA.5 sera", "XBB.2 sera")

{
  maps <- list(
    `discrete` = mapDisc,
    `continuous_fixtop_fixbottom` = mapContFixTopFixBottom,
    `continuous_fixtop` = mapContFixTop,
    `continuous_fixbottom` = mapContFixBottom,
    `continuous` = mapCont
  )
}

map_info <- long_map_list_info(maps) %>% filter(titer != "*")

homologous_ags <- c(
  "Beta sera discrete" = "Beta",
  "Delta sera discrete" = "Delta",
  "D614G sera discrete" = "D614G",
  "Alpha sera discrete" = "Alpha",
  "BA.1 sera discrete" = "BA.1",
  "BA.2 sera discrete" = "BA.2",
  "BA.2-12 sera discrete" = "BA.2-12",
  "B.1+E484K sera discrete" = "B.1+E484K",
  "BA.5 sera discrete" = "BA.5",
  "XBB.2 sera discrete" = "XBB.2",
  "Beta sera continuous_fixtop_fixbottom" = "Beta",
  "Delta sera continuous_fixtop_fixbottom" = "Delta",
  "D614G sera continuous_fixtop_fixbottom" = "D614G",
  "Alpha sera continuous_fixtop_fixbottom" = "Alpha",
  "BA.1 sera continuous_fixtop_fixbottom" = "BA.1",
  "BA.2 sera continuous_fixtop_fixbottom" = "BA.2",
  "BA.2-12 sera continuous_fixtop_fixbottom" = "BA.2-12",
  "B.1+E484K sera continuous_fixtop_fixbottom" = "B.1+E484K",
  "BA.5 sera continuous_fixtop_fixbottom" = "BA.5",
  "XBB.2 sera continuous_fixtop_fixbottom" = "XBB.2",
  "Beta sera continuous_fixtop" = "Beta",
  "Delta sera continuous_fixtop" = "Delta",
  "D614G sera continuous_fixtop" = "D614G",
  "Alpha sera continuous_fixtop" = "Alpha",
  "BA.1 sera continuous_fixtop" = "BA.1",
  "BA.2 sera continuous_fixtop" = "BA.2",
  "BA.2-12 sera continuous_fixtop" = "BA.2-12",
  "B.1+E484K sera continuous_fixtop" = "B.1+E484K",
  "BA.5 sera continuous_fixtop" = "BA.5",
  "XBB.2 sera continuous_fixtop" = "XBB.2",
  "Beta sera continuous_fixbottom" = "Beta",
  "Delta sera continuous_fixbottom" = "Delta",
  "D614G sera continuous_fixbottom" = "D614G",
  "Alpha sera continuous_fixbottom" = "Alpha",
  "BA.1 sera continuous_fixbottom" = "BA.1",
  "BA.2 sera continuous_fixbottom" = "BA.2",
  "BA.2-12 sera continuous_fixbottom" = "BA.2-12",
  "B.1+E484K sera continuous_fixbottom" = "B.1+E484K",
  "BA.5 sera continuous_fixbottom" = "BA.5",
  "XBB.2 sera continuous_fixbottom" = "XBB.2",
  "Beta sera continuous" = "Beta",
  "Delta sera continuous" = "Delta",
  "D614G sera continuous" = "D614G",
  "Alpha sera continuous" = "Alpha",
  "BA.1 sera continuous" = "BA.1",
  "BA.2 sera continuous" = "BA.2",
  "BA.2-12 sera continuous" = "BA.2-12",
  "B.1+E484K sera continuous" = "B.1+E484K",
  "BA.5 sera continuous" = "BA.5",
  "XBB.2 sera continuous" = "XBB.2"
)


### Make fold drop plots from real fold changes

# Compute fold drops for map
foldchange_table <- calculate_fold_change(map_info, sr_groups, homologous_ags)

plots <- list()
for(sr_g in sr_groups) {
  subset(foldchange_table, sr_group == sr_g & ag_name != "BA.2-2") %>%
    ggplot(
      aes(
        x = ag_name,
        y = mean_diff,
        ymin = mean_diff_lower,
        ymax = mean_diff_upper,
        color = map
      )
    ) + 
    geom_pointrange(
      position = position_dodge2(width = 0.7),
      size = 0.5
    ) +
    scale_y_continuous(
      breaks = function(x){
        ceiling(x[1]):floor(x[2])
      }
    ) +
    coord_cartesian(
      ylim = c(-8, 2)
    ) +
    scale_color_manual(
      values = list(
        `continuous_fixtop_fixbottom` = "firebrick4",
        `continuous_fixtop` = "red",
        `continuous_fixbottom` = "orange",
        `continuous` = "lightcoral",
        `discrete` = "blue"
      )
    ) +
    geom_hline(yintercept = 0,
               color = "#333333",
               linewidth = 0.1,
               linetype = "dashed"
    ) +
    labs(
      title = sr_g,
      x = "",
      y = "",
      color = "Method"
    ) +
    titerplot_theme() +
    scale_x_discrete(
      limits = rev(subset(foldchange_table, sr_group == sr_g & map == "discrete" & ag_name != "BA.2-2")$ag_name[order(subset(foldchange_table, sr_group == sr_g & map == "discrete" & ag_name != "BA.2-2")$mean_diff, na.last = NA)])
    ) -> gp
  plots <- c(plots, list(gp))
}

design <- "
  ABC
  DEF
  GHI
  J##
"

gps <- plots[[3]] + plots[[4]] +  plots[[2]] + plots[[1]]+ plots[[8]] + plots[[5]] + plots[[6]] + plots[[7]] + plots[[9]] + plots[[10]] + plot_layout(guides = "collect", design = design) + labs(y = "Fold change")
plot(gps)

ggsave("figures/fig_s4_discrete_and_continuous_prnt90_fold_drop/discrete_and_continuous_prnt90_fold_drop_update_240509.png", gps, width = 15, height = 11.5)


