
rm(list = ls())

library(tidyverse)
library(patchwork)
library(Racmacs)

source("Rcivaclib/map_longinfo.R")
source("Rcivaclib/diagnostics.R")
source("Rcivaclib/imputation.R")
source("Rcivaclib/scales.R")
source("Rcivaclib/common.R")

map <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

# Set serum group order
srGroups(map) <- factor(srGroups(map), sr_group_order)

residual_data <- residualErrorTable(map)

# Add antigen and sera info
residual_data %>% 
  left_join(long_ag_info(map), by = "ag_num") %>%
  left_join(long_sr_info(map), by = "sr_num") -> residual_data

source("Rcivaclib/diagnostics.R")
# Impute data
residual_data %>% 
  filter(
    measured_titer != "*"
  ) %>%
  filter(
    !paste(sr_group, ag_name) %in% c("BA.2 sera BQ.1.18", "BA.2 sera D614G",
                                     "BA.2 sera B.1+E484K", "BA.1 sera BQ.1.18",
                                     "BA.1 sera D614G", "BA.2-12 sera D614G",
                                     "D614G sera BF.7", "D614G sera BQ.1.18",
                                     "D614G sera XBB.2", "BA.5 sera BA.1",
                                     "BA.5 sera D614G", "XBB.2 sera D614G")
  ) %>%
  group_by(
    ag_name,
    sr_group
  ) %>%
  group_modify(
    .f = impute_table_residuals
  ) %>%
  filter(
    residual_error_imputed != Inf
  ) -> residual_data


# Plot the residual error with imputed values for < cases
## Boxplots split by antigens and sera
residual_data %>%
  mutate(
    ag_name = factor(ag_name, levels = c("D614G", "Alpha", "Delta", "Beta", "Mu", "B.1+E484K", "BA.2-12", "BA.2", "BA.1", "BA.4", "BA.5", "BF.7", "BQ.1.18", "BN.1.3.1", "XBB.2", "EG.5.1", "JN.1"))
  ) %>%
  ggplot(
    aes(
      x = ag_name,
      y = residual_error_imputed,
      fill = ag_name
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_boxplot(
    lwd = 0.25,
    outlier.shape = NA,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  ggbeeswarm::geom_beeswarm(
    aes(
      color = titer_type
    ),
    size = 0.5,
    cex = 0.5,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  facet_wrap(
    vars(sr_group)
  ) +
  scale_color_manual(
    values = c(
      "=" = "grey20",
      "<" = "lightblue",
      ">" = "pink"
    )
  ) +
  scale_fill_manual(
    values = agFillScale(map)
  ) +
  # scale_x_discrete(
  #   limits = agNames(map)
  # ) +
  coord_cartesian(
    ylim = c(-5, 5)
  ) +
  titerplot_theme() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(
    x = "",
    y = "Measured log titer - Fitted log titer"
  ) -> gp

ggsave("figures/fig_s11-12_map_vs_table_distances/map_vs_table_distances_by_sr_ag_240123.png", plot = gp, width = 8, height=6)

# Calculate summary parameters for the fit
residual_data %>% 
  ungroup() %>%
  filter(
    measured_titer != "*"
  ) %>%
  summarise(
    n = sum(measured_titer != "*"),
    fit = list(fit_censored_normal(
      lower_lims = residual_lower,
      upper_lims = residual_upper,
      mean = mean(residual_upper, na.rm = T),
      sd = sd(residual_upper, na.rm = T)
    )),
    .groups = "keep"
  ) %>% 
  mutate(
    mean = vapply(fit, \(x) x$estimate["mean"], numeric(1)),
    sd   = vapply(fit, \(x) x$estimate["sd"], numeric(1))
  ) -> residual_summary

fit_0 <- fit_censored_normal(
  lower_lims = residual_data$residual_lower,
  upper_lims = residual_data$residual_upper,
  mean = mean(0, na.rm = T),
  sd = sd(residual_data$residual_error_imputed, na.rm = T)
)

sprintf("Standard deviation with mean of 0: %s", fit_0$estimate["sd"])
sprintf("Standard deviation for error in both titers: %s", sqrt(fit_0$estimate["sd"] ^ 2 / 2))


# Plot a histogram with the overall distribution of the residuals
residual_data %>% 
  ggplot() + 
  geom_histogram(
    aes(
      x = residual_error_imputed,
      fill = titer_type
    ),
    binwidth = 0.2,
    color = "grey50",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = residual_summary$mean,
    lty = 2,
    lwd = 0.6
  ) +
  coord_cartesian(
    xlim = c(-8, 8)
  ) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  labs(
    x = "Measured log titer - Fitted log titer",
    y = "Count",
    fill = "Titer type",
    subtitle = sprintf(
      "µ = %s, σ = %s",
      round(residual_summary$mean, 2),
      round(residual_summary$sd, 2)
    )
  ) -> gp_hist


# Plot the fitted against the measured titers as a scatter plot
ggplot2::ggplot(
  data = residual_data) +
  ggplot2::geom_point(
    mapping = ggplot2::aes(
      y     = measured_logtiter_upper,
      x     = predicted_logtiter,
      color = titer_type
    ),
    alpha = 0.4
  ) +
  ggplot2::theme_bw() +
  coord_cartesian(
    xlim = c(1, 9),
    ylim = c(1, 9)
  ) +
  scale_x_continuous(
    breaks = 1:9,
    labels = c("20", 2^(2:9)*10)
  ) +
  scale_y_continuous(
    breaks = 1:9,
    labels = c("20", 2^(2:9)*10)
  ) +
  ylab("Measured titers") +
  xlab("Fitted titers") +
  geom_abline(intercept = 0, lty=2, col = "grey") +
  guides(color = "none") +
  theme(
    text = element_text(size=15),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) -> gp_scatter

gp_scatter

plots <- gp_scatter + gp_hist + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 25))

plots

ggsave(plot = plots, "figures/fig_s11-12_map_vs_table_distances/map_vs_table_distances_scatter_hist_240123.png", width=10, height=5)


sqrt(sum(residual_data$residual_error_imputed^2) / length(residual_data$residual_error_imputed))


# Calculate statistics for detectable titers only
residual_data %>% 
  ungroup() %>%
  filter(
    titer_type == "="
  ) %>%
  summarise(
    n = sum(titer_type == "="),
    fit = list(fit_censored_normal(
      lower_lims = residual_lower,
      upper_lims = residual_upper,
      mean = mean(residual_upper, na.rm = T),
      sd = sd(residual_upper, na.rm = T)
    )),
    .groups = "keep"
  ) %>% 
  mutate(
    mean = vapply(fit, \(x) x$estimate["mean"], numeric(1)),
    sd   = vapply(fit, \(x) x$estimate["sd"], numeric(1))
  ) -> residual_summary_detectable

residual_data %>% 
  ungroup() %>%
  filter(
    titer_type == "="
  ) -> residual_data_detectable

fit_0_detectable <- fit_censored_normal(
  lower_lims = residual_data_detectable$residual_lower,
  upper_lims = residual_data_detectable$residual_upper,
  mean = mean(0, na.rm = T),
  sd = sd(residual_data_detectable$residual_error_imputed, na.rm = T)
)

sprintf("Standard deviation with mean of 0: %s", fit_0_detectable$estimate["sd"])
sprintf("Standard deviation for error in both titers: %s", sqrt(fit_0_detectable$estimate["sd"] ^ 2 / 2))

