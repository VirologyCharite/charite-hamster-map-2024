
rm(list = ls())

library(tidyverse)
library(Racmacs)
library(patchwork)
library(titertools)

source("Rcivaclib/scales.R")
source("Rcivaclib/calculate_fold_change.R")
source("Rcivaclib/map_longinfo.R")

map50 <- read.acmap("data/240123-hamster/map-discrete-50.ace")
map75 <- read.acmap("data/240123-hamster/map-discrete-75.ace")
map90 <- read.acmap("data/240123-hamster/map-discrete-90.ace")
map99 <- read.acmap("data/240123-hamster/map-discrete-99.ace")
map50cont <- read.acmap("data/240123-hamster/map-continuous-fixtop-fixbottom-50.ace")
map75cont <- read.acmap("data/240123-hamster/map-continuous-fixtop-fixbottom-75.ace")
map90cont <- read.acmap("data/240123-hamster/map-continuous-fixtop-fixbottom-90.ace")
map99cont <- read.acmap("data/240123-hamster/map-continuous-fixtop-fixbottom-99.ace")

sr_groups <- c("Beta sera", "Delta sera", "D614G sera", "Alpha sera",
               "BA.1 sera", "BA.2 sera", "BA.2-12 sera", "B.1+E484K sera",
               "BA.5 sera", "XBB.2 sera")

{
  maps <- list(
    `discrete_nt50` = map50,
    `discrete_nt75` = map75,
    `discrete_nt90` = map90,
    `discrete_nt99` = map99,
    `continuous_nt50` = map50cont,
    `continuous_nt75` = map50cont,
    `continuous_nt90` = map90cont,
    `continuous_nt99` = map90cont
  )
}


map_info <- long_map_list_info(maps) %>% filter(titer != "*")

homologous_ags <- c(
  "Beta sera discrete_nt50" = "Beta",
  "Delta sera discrete_nt50" = "Delta",
  "D614G sera discrete_nt50" = "D614G",
  "Alpha sera discrete_nt50" = "Alpha",
  "BA.1 sera discrete_nt50" = "BA.1",
  "BA.2 sera discrete_nt50" = "BA.2",
  "BA.2-12 sera discrete_nt50" = "BA.2-12",
  "B.1+E484K sera discrete_nt50" = "B.1+E484K",
  "BA.5 sera discrete_nt50" = "BA.5",
  "XBB.2 sera discrete_nt50" = "XBB.2",
  "Beta sera discrete_nt75" = "Beta",
  "Delta sera discrete_nt75" = "Delta",
  "D614G sera discrete_nt75" = "D614G",
  "Alpha sera discrete_nt75" = "Alpha",
  "BA.1 sera discrete_nt75" = "BA.1",
  "BA.2 sera discrete_nt75" = "BA.2",
  "BA.2-12 sera discrete_nt75" = "BA.2-12",
  "B.1+E484K sera discrete_nt75" = "B.1+E484K",
  "BA.5 sera discrete_nt75" = "BA.5",
  "XBB.2 sera discrete_nt75" = "XBB.2",
  "Beta sera discrete_nt90" = "Beta",
  "Delta sera discrete_nt90" = "Delta",
  "D614G sera discrete_nt90" = "D614G",
  "Alpha sera discrete_nt90" = "Alpha",
  "BA.1 sera discrete_nt90" = "BA.1",
  "BA.2 sera discrete_nt90" = "BA.2",
  "BA.2-12 sera discrete_nt90" = "BA.2-12",
  "B.1+E484K sera discrete_nt90" = "B.1+E484K",
  "BA.5 sera discrete_nt90" = "BA.5",
  "XBB.2 sera discrete_nt90" = "XBB.2",
  "Beta sera discrete_nt99" = "Beta",
  "Delta sera discrete_nt99" = "Delta",
  "D614G sera discrete_nt99" = "D614G",
  "Alpha sera discrete_nt99" = "Alpha",
  "BA.1 sera discrete_nt99" = "BA.1",
  "BA.2 sera discrete_nt99" = "BA.2",
  "BA.2-12 sera discrete_nt99" = "BA.2-12",
  "B.1+E484K sera discrete_nt99" = "B.1+E484K",
  "BA.5 sera discrete_nt99" = "BA.5",
  "XBB.2 sera discrete_nt99" = "XBB.2",
  "Beta sera continuous_nt50" = "Beta",
  "Delta sera continuous_nt50" = "Delta",
  "D614G sera continuous_nt50" = "D614G",
  "Alpha sera continuous_nt50" = "Alpha",
  "BA.1 sera continuous_nt50" = "BA.1",
  "BA.2 sera continuous_nt50" = "BA.2",
  "BA.2-12 sera continuous_nt50" = "BA.2-12",
  "B.1+E484K sera continuous_nt50" = "B.1+E484K",
  "BA.5 sera continuous_nt50" = "BA.5",
  "XBB.2 sera continuous_nt50" = "XBB.2",
  "Beta sera continuous_nt90" = "Beta",
  "Delta sera continuous_nt90" = "Delta",
  "D614G sera continuous_nt90" = "D614G",
  "Alpha sera continuous_nt90" = "Alpha",
  "BA.1 sera continuous_nt90" = "BA.1",
  "BA.2 sera continuous_nt90" = "BA.2",
  "BA.2-12 sera continuous_nt90" = "BA.2-12",
  "B.1+E484K sera continuous_nt90" = "B.1+E484K",
  "BA.5 sera continuous_nt90" = "BA.5",
  "XBB.2 sera continuous_nt90" = "XBB.2",
  "Beta sera continuous_nt75" = "Beta",
  "Delta sera continuous_nt75" = "Delta",
  "D614G sera continuous_nt75" = "D614G",
  "Alpha sera continuous_nt75" = "Alpha",
  "BA.1 sera continuous_nt75" = "BA.1",
  "BA.2 sera continuous_nt75" = "BA.2",
  "BA.2-12 sera continuous_nt75" = "BA.2-12",
  "B.1+E484K sera continuous_nt75" = "B.1+E484K",
  "XBB.2 sera continuous_nt75" = "XBB.2",
  "BA.5 sera continuous_nt75" = "BA.5",
  "Beta sera continuous_nt99" = "Beta",
  "Delta sera continuous_nt99" = "Delta",
  "D614G sera continuous_nt99" = "D614G",
  "Alpha sera continuous_nt99" = "Alpha",
  "BA.1 sera continuous_nt99" = "BA.1",
  "BA.2 sera continuous_nt99" = "BA.2",
  "BA.2-12 sera continuous_nt99" = "BA.2-12",
  "B.1+E484K sera continuous_nt99" = "B.1+E484K",
  "BA.5 sera continuous_nt99" = "BA.5",
  "XBB.2 sera continuous_nt99" = "XBB.2"
)


### Make fold drop plots from real fold changes

# Compute fold drops for map
foldchange_table <- calculate_fold_change(map_info, sr_groups, homologous_ags)


# Plots for discrete titers
plots <- list()
for(sr_g in sr_groups) {
  subset(foldchange_table, sr_group == sr_g & ag_name != "BA.2-2") %>%
    filter(map %in% c("discrete_nt50", "discrete_nt75", "discrete_nt90", "discrete_nt99")) %>%
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
        `discrete_nt50` = "firebrick4",
        `discrete_nt75` = "red",
        `discrete_nt90` = "orange",
        `discrete_nt99` = "lightcoral"
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
      color = "Sensitivity level"
    ) +
    titerplot_theme() +
    scale_x_discrete(
      limits = rev(subset(foldchange_table, sr_group == sr_g & map == "discrete_nt90" & ag_name != 'BA.2-2')$ag_name[order(subset(foldchange_table, sr_group == sr_g & map == 'discrete_nt90' & ag_name != 'BA.2-2')$mean_diff, na.last = NA)])
    ) -> gp
  plots <- c(plots, list(gp))
}

design <- "
  ABC
  DEF
  GHI
  J##
"

gps <- plots[[3]] + plots[[4]] +  plots[[2]] + plots[[1]]+ plots[[8]] + plots[[5]] + plots[[6]] + plots[[7]] + plots[[9]] + plots[[10]] + plot_layout(guides = "collect", design = design) + labs(y = 'Fold change')
plot(gps)

ggsave("figures/fig_s1-2_sensitivity_levels_fold_drop/sensitivity_levels_fold_drop_discrete_updated_240509.png", gps, width = 15, height = 11.5)



# Plots for continuous titers
plots <- list()
for(sr_g in sr_groups) {
  subset(foldchange_table, sr_group == sr_g & ag_name != "BA.2-2") %>%
    filter(map %in% c("continuous_nt50", "continuous_nt75", "continuous_nt90", "continuous_nt99")) %>%
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
        `continuous_nt50` = "firebrick4",
        `continuous_nt75` = "red",
        `continuous_nt90` = "orange",
        `continuous_nt99` = "lightcoral"
      )
    ) +
    geom_hline(yintercept = 0,
               color = "#333333",
               size = 0.1,
               linetype = "dashed"
    ) +
    labs(
      title = sr_g,
      x = "",
      y = "",
      color = "Sensitivity level"
    ) +
    titerplot_theme() +
    scale_x_discrete(
      limits = rev(subset(foldchange_table, sr_group == sr_g & map == "continuous_nt90" & ag_name != "BA.2-2")$ag_name[order(subset(foldchange_table, sr_group == sr_g & map == "continuous_nt90" & ag_name != "BA.2-2")$mean_diff, na.last = NA)])
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

ggsave("figures/sensitivity_levels_fold_drop/sensitivity_levels_fold_drop_continuous_updated_240509.png", gps, width = 15, height = 11.5)
