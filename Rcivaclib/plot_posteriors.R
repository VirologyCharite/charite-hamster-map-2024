# Plot a bar plot of the samples from the posterior distributions
# data: The calculated highest posterior density intervals.
# map_order: The order in which the datasets should be plotted.
# map_cols: The colors in which the datasets should be plotted.
# xlabel: The label for the x-axis label.
plot_posteriors_hpdi <- function(data, map_order, map_cols, xlabel) {

  data %>%
    ggplot(
      aes(
        y = map,
        x = mean_,
        xmin = ci_low,
        xmax = ci_high,
        color = map
      )
    ) +
    geom_linerange(
      size = 5,
      alpha = 0.5
    ) +
    geom_point(
      size = 4,
      alpha = 0.1
    ) +
    geom_point(
      size = 3,
      color = 'white'
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 'dashed',
      color = 'black',
      alpha = 0.5
    ) +
    scale_y_discrete(
      limits = map_order
  ) +
  labs(
    x = xlabel,
    y = ''
  ) +
  scale_color_manual(
    values = map_cols
  ) +
  coord_cartesian(
    xlim = c(-10, 10)
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 15)) +
  guides(color = 'none', size = 'none') -> gp
}
