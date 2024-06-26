
library(Racmacs)
library(tidyverse)
library(dplyr) 
library(patchwork)

library(titertools)

source("Rcivaclib/map_longinfo.R")
source("Rcivaclib/gmts.R")
source("Rcivaclib/plot_posteriors.R")

sr_groups <- c("D614G sera", "Alpha sera", "Beta sera")

# Load maps
hamster_subset <- read.acmap("figures/fig_3_hamster_human_comparison/hamster_subset.ace")
srOutlineWidth(hamster_subset) <- 1
human_subset <- read.acmap("figures/fig_3_hamster_human_comparison/human_subset.ace")

hamsterP <- procrustesMap(hamster_subset, human_subset)
humanP <- procrustesMap(human_subset, hamster_subset)

## Plot titers
data_gmts <- readRDS("figures/fig_3_hamster_human_comparison/titer_magnitude/gmts.rds")

data_gmts %>%
  mutate(
    sr_group = factor(sr_group, levels = c("D614G sera", "Alpha sera", "Beta sera"))
  ) %>%
  ggplot(
    aes(
      x = ag_name,
      y = gmt,
      ymin = gmt_upper,
      ymax = gmt_lower,
      color = map,
      group = map
    )
  ) + 
  geom_pointrange(
    position = position_dodge2(
      width = 0.8
    )
  ) +
  geom_line(
    linetype = "solid", position = position_dodge2(
      width = 0.8
    )
  ) +
  scale_y_continuous(breaks=seq(0, 10),
                     labels=c("<20", "20", "40", "80", "160", "320", "640", "1280", "2560", "5120", "10240"),
                     minor_breaks=F) +
  scale_color_manual(
    values = c("#235c3f", "#d95f02"),
    name = "Dataset"
  ) +
  theme_bw() +
  labs(
    x = "",
    y = "Titer"
  ) +
  theme(
    text = element_text(size = 15),
    #legend.position = "top",
    legend.key = element_rect(fill = "white"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.spacing.x = unit(0.85, "lines")
  ) +
  guides(color = "none") +
  facet_wrap(~sr_group) -> gp_titer_magnitude


## Do the plotting of estimated titer magnitude differences
# Read in the data
draws <- readRDS("figures/fig_3_hamster_human_comparison/titer_magnitude/dataset_magnitude_effect_posterior_samples.rds")

# Adjust the column headers
colnames(draws) <- c("hamster", "human", ".chain", ".iteration", ".draw")

draws$differences <- draws$human - draws$hamster

# Convert to long format
draws %>%
  pivot_longer(
    cols = c("hamster", "human", "differences"),
    names_to = "map") -> draws_long

# Calculate the 95% HPDI
par_ci <- bayestestR::ci(draws, ci = 0.95, method = "HDI")

# Convert the calculated 95% HPDI to a dataframe for plotting
data <- tibble(map = par_ci$Parameter, ci_low = par_ci$CI_low, ci_high = par_ci$CI_high, mean_ = colMeans(draws)[ -c(3:5) ])

plot_posteriors_hpdi(data, c("hamster", "human", "differences"), c(`hamster` = "#235c3f", `human` = "#d95f02", `differences` = "grey30"),
                     xlabel = "") +
  coord_cartesian(xlim = c(-10, 10)) +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size=12, hjust = 0.5, margin=margin(0,0,-20,0)),
    axis.title.x = element_text(size = 12, margin = margin(t = -50, r = 0, b = 0, l = 0))
  ) +
  labs(title  = "Estimated magnitude\neffect", x = "Magnitude effect") +
  scale_y_discrete(
    limits = rev(c("hamster", "human", "differences")),
    labels = rev(c("Hamster", "Human", "Human - \n Hamster"))
  ) -> gp_dataset_effect
gp_dataset_effect


## Plot the fold change
foldchange_table <- readRDS("figures/fig_3_hamster_human_comparison/foldchange/foldchange_table.rds")

plots <- list()
for(sr_g in sr_groups) {
  subset(foldchange_table, sr_group == sr_g) %>%
    ggplot(
      aes(
        x = ag_name,
        y = mean_diff,
        ymin = mean_diff_lower,
        ymax = mean_diff_upper,
        color = map,
        group = map
      )
    ) + 
    geom_hline(yintercept = 0,
               color = "#333333",
               size = 0.1
    ) +
    geom_pointrange(
      position = position_dodge2(width = 0.5)
    ) +
    geom_line(
      position = position_dodge2(width = 0.5),
      linetype = "solid"
    ) +
    scale_y_continuous(
      breaks = function(x){
        ceiling(x[1]):floor(x[2])
      }
    ) +
    coord_cartesian(
      ylim = c(-7, 1)
    ) +
    labs(
      title = sr_g,
      x = ",
      y = "
    ) +
    scale_color_manual(
      values = c("#235c3f", "#d95f02"),
      name = "Dataset"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(size=12, hjust = 0.5),
      legend.key = element_rect(fill = "white"),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    ) +
    guides(color = "none") +
    scale_x_discrete(
      limits = rev(subset(foldchange_table, sr_group == sr_g & map == "hamster")$ag_name[order(subset(foldchange_table, sr_group == sr_g & map == "hamster")$mean_diff, na.last = NA)])
    ) -> gp
  plots <- c(plots, list(gp))
}

foldchange_magnitude <- plots[[1]] + labs(y = "Foldchange") +
  plots[[2]] + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank()) +
  plots[[3]] + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank()) +
  plot_layout(guides = "collect")
plot(foldchange_magnitude)


## Plot foldchange differences
# Read in the data
draws <- readRDS("figures/fig_3_hamster_human_comparison/foldchange/slope_calculation_dataset_slope_effect_draws.rds")

# Adjust the column headers
colnames(draws) <- c("hamster", "human", ".chain", ".iteration", ".draw")

draws$differences <- draws$human - draws$hamster

# Convert to long format
draws %>%
  pivot_longer(
    cols = c("hamster", "human", "differences"),
    names_to = "map") -> draws_long

# Calculate the 95% HPDI
par_ci <- bayestestR::ci(draws, ci = 0.95, method = "HDI")

# Convert the calculated 95% HPDI to a dataframe for plotting
data <- tibble(map = par_ci$Parameter, ci_low = par_ci$CI_low, ci_high = par_ci$CI_high, mean_ = colMeans(draws)[ -c(3:5) ])

# Do the plotting
plot_posteriors_hpdi(data, c("hamster", "human", "differences"), c(`hamster` = "#235c3f", `human` = "#d95f02", `differences` = "grey30"),
                            xlabel = "") +
  coord_cartesian(xlim = c(0, 1.25)) +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size=12, hjust = 0.5),
    axis.title.x = element_text(size = 12, margin = margin(t = -50, r = 0, b = 0, l = 0))
  ) +
  labs(title  = "Estimated slope effect", x = "Slope effect") +
  scale_y_discrete(
    limits = rev(c("hamster", "human", "differences")),
    labels = rev(c("Hamster", "Human", "Human - \n Hamster"))
  ) -> gp_slope_diff
gp_slope_diff


hamsterMap <- ggplot(hamsterP, xlim = c(-3, 1), ylim = c(-2, 2), cex = 0.3) + 
  labs(title = "Hamster") +
  theme(
    plot.title = element_text(size=12, hjust = 0.5, margin=margin(0,0,-20,0))
  )
hamsterMap

humanMap <- ggplot(humanP, xlim = c(-3, 1), ylim = c(-2, 2), cex = 0.3) +
  labs(title = "Human") +
  theme(
    plot.title = element_text(size=12, hjust = 0.5, margin=margin(0,0,5,0))
  )
humanMap

# Assemble the figures
design <- "
  ABC
  DEF
"
xx <- gp_titer_magnitude + gp_dataset_effect + hamsterMap + foldchange_magnitude + gp_slope_diff + humanMap + 
  plot_layout(design = design, widths = c(3, 1, 2))
xx
ggsave("figures/fig_3_hamster_human_comparison/figure.pdf", width = 12, height = 8)

