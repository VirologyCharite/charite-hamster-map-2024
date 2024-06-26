
rm(list = ls())

library(tidyverse)
library(Racmacs)
library(patchwork)
library(titertools)

source("Rcivaclib/scales.R")
source("Rcivaclib/calculate_fold_change.R")
source("Rcivaclib/map_longinfo.R")

set.seed(100)

pre_omi_sr <- c("D614G sera", "Alpha sera", "Beta sera", "Delta sera",
                "B.1+E484K sera")

pre_omi_sr_no_484 <- c("D614G sera", "Alpha sera", "Beta sera", "Delta sera")

map <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

sr_groups <- c("Beta sera", "Delta sera", "D614G sera", "Alpha sera",
               "BA.1 sera", "BA.2 sera", "BA.2-12 sera", "B.1+E484K sera",
               "BA.5 sera", "XBB.2 sera")

agNames(map) <- ag_names
srGroups(map) <- sr_groups_ind

{
  maps <- list(
    `civ_prnt90` = map
  )
}

data <- long_map_list_info(maps)

homologous_ags <- c(
  "Beta sera civ_prnt90" = "Beta",
  "Delta sera civ_prnt90" = "Delta",
  "D614G sera civ_prnt90" = "D614G",
  "Alpha sera civ_prnt90" = "Alpha",
  "BA.1 sera civ_prnt90" = "BA.1",
  "BA.2 sera civ_prnt90" = "BA.2",
  "BA.2-12 sera civ_prnt90" = "BA.2-12",
  "B.1+E484K sera civ_prnt90" = "B.1+E484K",
  "BA.5 sera civ_prnt90" = "BA.5",
  "XBB.2 sera civ_prnt90" = "XBB.2"
)

foldchange_table <- calculate_fold_change(data, sr_groups, homologous_ags)
saveRDS(foldchange_table, "figures/fig_s8_fold_change/foldchange_table_240123.rds")
stop()

foldchange_table <- readRDS("figures/fig_s8_fold_change/foldchange_table_240123.rds")

# Convert log2 to fold change
foldchange <- \(x) {
  
  xabs <- abs(x)
  foldchange <- 2^xabs
  foldchange[x < 0] <- -foldchange[x < 0]
  as.character(round(foldchange, 1))
}

x <- tibble(ag_name = foldchange_table$ag_name, sr_group = foldchange_table$sr_group, adjusted_diff=sapply(foldchange_table$mean_diff, foldchange))


# Look up fold changes for results section
foldchange(min(subset(foldchange_table, ag_name == "BA.1" & sr_group %in% pre_omi_sr_no_484)$mean_diff))
foldchange(max(subset(foldchange_table, ag_name == "BA.1" & sr_group %in% pre_omi_sr_no_484)$mean_diff))

foldchange(min(subset(foldchange_table, ag_name == "BA.4" & sr_group %in% pre_omi_sr_no_484)$mean_diff))
foldchange(max(subset(foldchange_table, ag_name == "BA.4" & sr_group %in% pre_omi_sr_no_484)$mean_diff))

foldchange(min(subset(foldchange_table, ag_name == "BA.5" & sr_group %in% pre_omi_sr_no_484)$mean_diff))
foldchange(max(subset(foldchange_table, ag_name == "BA.5" & sr_group %in% pre_omi_sr_no_484)$mean_diff))

r <- titertools::log2diff(
  titers2 = subset(data, sr_group == "B.1+E484K convalescent" & ag_name == "BA.4")$titer,
  titers1 = subset(data, sr_group == "B.1+E484K convalescent" & ag_name == "BA.1")$titer,
  dilution_stepsize = 0,
  ci_method = "HDI",
  mu_prior_mu = 0,
  mu_prior_sigma = 100,
  sigma_prior_alpha = 2,
  sigma_prior_beta = 0.75
)
foldchange(r["mean", "estimate"])
foldchange(r["mean", "lower"])
foldchange(r["mean", "upper"])

titertools::log2diff(
  titers2 = subset(data, sr_group == "B.1+E484K convalescent" & ag_name == "BA.5")$titer,
  titers1 = subset(data, sr_group == "B.1+E484K convalescent" & ag_name == "BA.1")$titer,
  dilution_stepsize = 0,
  ci_method = "HDI",
  mu_prior_mu = 0,
  mu_prior_sigma = 100,
  sigma_prior_alpha = 2,
  sigma_prior_beta = 0.75
)

foldchange(subset(foldchange_table, ag_name == "Mu" & sr_group == "B.1+E484K convalescent")$mean_diff)
foldchange(subset(foldchange_table, ag_name == "Mu" & sr_group == "B.1+E484K convalescent")$mean_diff_lower)
foldchange(subset(foldchange_table, ag_name == "Mu" & sr_group == "B.1+E484K convalescent")$mean_diff_upper)


for (ag in unique(data$ag_name)) {
  print(ag)
  result <- titertools::log2diff(
    titers2 = subset(data, sr_group == "BA.1 convalescent" & ag_name == ag)$titer,
    titers1 = subset(data, sr_group == "BA.1 convalescent" & ag_name == "BA.5")$titer,
    dilution_stepsize = 0,
    ci_method = "HDI",
    mu_prior_mu = 0,
    mu_prior_sigma = 100,
    sigma_prior_alpha = 2,
    sigma_prior_beta = 0.75
  )
  print(result)
}


for (ag in unique(data$ag_name)) {
  print(ag)
  result <- titertools::log2diff(
    titers2 = subset(data, sr_group == "BA.1 convalescent" & ag_name == ag)$titer,
    titers1 = subset(data, sr_group == "BA.1 convalescent" & ag_name == "BA.4")$titer,
    dilution_stepsize = 0,
    ci_method = "HDI",
    mu_prior_mu = 0,
    mu_prior_sigma = 100,
    sigma_prior_alpha = 2,
    sigma_prior_beta = 0.75
  )
  print(result)
}

titertools::log2diff(
  titers2 = subset(data, sr_group %in% c("BA.2 convalescent", "BA.2 convalescent-12") & ag_name == "BA.1")$titer,
  titers1 = subset(data, sr_group %in% c("BA.2 convalescent", "BA.2 convalescent-12") & ag_name == "BA.2")$titer,
  dilution_stepsize = 0,
  ci_method = "HDI",
  mu_prior_mu = 0,
  mu_prior_sigma = 100,
  sigma_prior_alpha = 2,
  sigma_prior_beta = 0.75
)

titertools::log2diff(
  titers2 = subset(data, sr_group %in% c("BA.2 convalescent", "BA.2 convalescent-12") & ag_name == "BA.4")$titer,
  titers1 = subset(data, sr_group %in% c("BA.2 convalescent", "BA.2 convalescent-12") & ag_name == "BA.2")$titer,
  dilution_stepsize = 0,
  ci_method = "HDI",
  mu_prior_mu = 0,
  mu_prior_sigma = 100,
  sigma_prior_alpha = 2,
  sigma_prior_beta = 0.75
)

titertools::log2diff(
  titers2 = subset(data, sr_group %in% c("BA.2 convalescent", "BA.2 convalescent-12") & ag_name == "BA.5")$titer,
  titers1 = subset(data, sr_group %in% c("BA.2 convalescent", "BA.2 convalescent-12") & ag_name == "BA.2")$titer,
  dilution_stepsize = 0,
  ci_method = "HDI",
  mu_prior_mu = 0,
  mu_prior_sigma = 100,
  sigma_prior_alpha = 2,
  sigma_prior_beta = 0.75
)

# Plot fold change

colorListNew <- agFill(map)
names(colorListNew) <- agNames(map)

plots <- list()
for(sr_g in sr_groups) {
  subset(foldchange_table, sr_group == sr_g & ag_name != "BA.2-2") %>%
    ggplot(
      aes(
        x = ag_name,
        y = mean_diff,
        ymin = mean_diff_lower,
        ymax = mean_diff_upper,
        #color = map
        color = ag_name
      )
    ) + 
    geom_pointrange(
      position = position_dodge2(width = 0.5),
      size = 0.8
    ) +
    scale_y_continuous(
      breaks = function(x){
        ceiling(x[1]):floor(x[2])
      }
    ) +
    coord_cartesian(
      ylim = c(-8, 2)
      #ylim = c(-11, 4)
    ) +
    scale_color_manual(
      values = colorListNew
    ) +
    geom_hline(yintercept = 0,
               color = "#333333",
               linewidth = 0.1,
               linetype = "dashed"
    ) +
    labs(
      title = sr_g,
      x = "",
      y = ""
    ) +
    titerplot_theme() +
    theme(
      plot.title = element_text(size=12),
      axis.text.x = element_text(size=10)
    ) +
    guides(color = "none") +
    scale_x_discrete(
      limits = rev(subset(foldchange_table, sr_group == sr_g & map == "civ_prnt90" & ag_name != "BA.2-2")$ag_name[order(subset(foldchange_table, sr_group == sr_g & map == "civ_prnt90" & ag_name != "BA.2-2")$mean_diff, na.last = F)])
    ) -> gp
  plots <- c(plots, list(gp))
}

design <- "
  ABCDE
  FGHIJ
"

gps <- plots[[3]] + labs(y = "Fold change") + plots[[4]] +  plots[[2]] + plots[[1]] + plots[[8]] + plots[[6]] + labs(y = "Fold change") + plots[[7]] + plots[[5]] + plots[[9]] + plots[[10]] + plot_layout(guides = "collect", design = design)

plot(gps)

ggsave("figures/fig_s8_fold_change/fold_change_240509.png", gps, width = 17, height = 6)


