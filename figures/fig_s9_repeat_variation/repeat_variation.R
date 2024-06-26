
# Setup workspace

rm(list = ls())
library(tidyverse)

data <- read.csv(file = "figures/fig_s9_repeat_variation/data_240123.csv")


data %>%
  ggplot(
    aes(
      x = prnt50ContFixbottom
    )
  ) +
  geom_histogram(
    bins = 100,
    fill = "darkblue",
    alpha = 0.5
  ) +
  xlim(
    c(-2.5, 2.5)
  ) +
  geom_vline(
    xintercept = mean(subset(data, prnt50ContFixbottom > -2.5)$prnt50ContFixbottom),
    linetype = "dashed"
  ) +
  labs(
    x = "log titer run 1 - log titer run 2",
    y = "Count",
    title = sprintf(
      "µ = %s, σ = %s",
      round(mean(subset(data, prnt50ContFixbottom > -2.5)$prnt50ContFixbottom), 2),
      round(sd(subset(data, prnt50ContFixbottom > -2.5)$prnt50ContFixbottom), 2)
    )
  ) +
  theme_bw()
ggsave("figures/fig_s9_repeat_variation/repeat_variation_240123.png", width = 5, height = 4)

# Standard deviation accounting for measurement error in the first and second repeat
sqrt(sum((subset(data, prnt50ContFixbottom > -2.5)$prnt50ContFixbottom)^2) / length(subset(data, prnt50ContFixbottom > -2.5)$prnt50ContFixbottom))

