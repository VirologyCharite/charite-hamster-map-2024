
rm(list = ls())

library(tidyverse)
library(Racmacs)
library(covutils)
library(patchwork)
library(titertools)

source("Rcivaclib/common.R")

map <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

agNames(map) <- c("BF.7", "BN.1.3.1", "BQ.1.18", "EG.5.1", "JN.1", "Alpha", "BA.1", "BA.2-12",
                  "BA.2", "BA.4", "BA.5", "Beta",
                  "Delta", "Mu", "D614G", "XBB.2", "B.1+E484K")

data <- long_map_info(map)

data$logtiter[data$titertype == 2] <- data$logtiter[data$titertype == 2] + 1
data$logtiter[data$titertype == 3] <- data$logtiter[data$titertype == 3] - 1

titles <- data.frame(
  row.names = c("BA.2-26729_2", "Alpha-21528", "Beta-22131", "Delta-25853_19", "BA.1-26335",
                "BA.2-26729_12", "B.1+E484K", "WT-984", "BA.5", "XBB.2"),
  val = c("BA.2 sera", "Alpha sera", "Beta sera", "Delta sera",
          "BA.1 sera", "BA.2-12 sera", "B.1+E484K sera", "D614G sera",
          "BA.5 sera", "XBB.2 sera")
)

# Add coloured boxes
add_coloured_boxes <- function(fig) {
  for (n in seq_len(numAntigens(map))) {
    fig <- fig +
      annotate(
        "tile",
        x = agNames(map)[n],
        y = 0,
        height = 1,
        fill = agFill(map)[n],
        color = NA,
        alpha = 0.9
      )
  }
  fig
}



do_plot <- function(d) {
  d %>%
    mutate(
      logtiter = ifelse(d$titertype == 2, NaN, d$logtiter)
    ) -> dnd
  d %>%
    ggplot(
      aes(
        x = ag_name,
        y = logtiter,
        color = sr_group,
        group = sr_name,
        shape = titertype,
        alpha = titertype
      )
    ) +
    geom_line(
      data = dnd,
      alpha = 1,
      size = 0.7
    ) +
    geom_line(
      alpha = 0.5,
      size = 0.2,
      linetype="42"
    ) +
    coord_cartesian(
      ylim = c(0, 9.5),
      xlim = c(1.1, 16.9)
    ) +
    geom_point(
      data = dnd,
      size = 2,
      alpha = 1
    ) +
    geom_point(
      size = 2,
      alpha = 0.2
    ) +
    theme_mapplot() +
    scale_y_continuous(
      breaks = seq(1, 9),
      labels = c("20", "40", "80", "160", "320", "640", "1280", "2560", "5120"),
      minor_breaks = NULL
    ) +
    scale_shape_manual(
      values = c(
        `1` = 16,
        `2` = 21,
        `3` = 21
      )
    ) +
    scale_color_discrete(
      name = "Serum"
    ) +
    guides(shape="none", color = "none", linetype = "none") +
    scale_color_manual(
      values = sr_group_cols
    ) +
    theme(
      plot.title = element_text(size=12)
    ) +
    labs(
      x = "",
      y = "",
      title = sg
    ) -> gp
  gp <- add_coloured_boxes(gp)
  
  gp
}


plots <- list()
for (sg in unique(data$sr_group)) {
  data %>%
    filter(sr_group == sg) %>%
    mutate(
      ag_name = factor(ag_name, levels = c("D614G", "Alpha", "Delta", "Beta", "Mu",
                                           "B.1+E484K", "BA.2", "BA.2-12", "BA.1", "BA.4", "BA.5", 
                                           "BF.7", "BQ.1.18", "BN.1.3.1", "XBB.2", "EG.5.1",
                                           "JN.1"))
    ) -> d
    do_plot(d) -> gp
  
  plots <- c(plots, list(gp))
}



wt <- plots[[9]] + theme(axis.text.x = element_blank(),
                          axis.ticks.x = element_blank())


b117 <- plots[[3]] + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank())

b16172 <- plots[[5]] + theme(axis.text.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.title.y = element_blank(),
                             axis.text.x = element_blank(),
                             axis.ticks.x = element_blank())

b1351 <- plots[[4]] + theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank())

e484k <- plots[[8]] + theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank())

ba2 <- plots[[1]] + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank() )

ba212 <- plots[[7]]

ba1 <- plots[[6]] + theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title.y = element_blank() )

ba5 <- plots[[10]] + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank() )

xbb2 <- plots[[2]] + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank())

design <- "
  ABCDE
  FGHIJ
"

gps <- wt + b117 + b16172 + b1351 + e484k + ba212 + ba2 + ba1 + ba5 + xbb2 + plot_layout(guides = "collect", design = design)

plot(gps)

ggsave("figures/fig_1_titerplot/titerplot_240509.pdf", gps, width = 15, height = 6)
ggsave("figures/fig_1_titerplot/titerplot_240509.png", gps, width = 15, height = 6)


# Print the GMTs mentioned in the text
get_titer <- function(titers, what = "mean") {
  logtiter <- titertools::gmt(
    titers = titers,
    ci_method = "HDI",
    sigma_prior_alpha = 2,
    sigma_prior_beta = 0.75,
    mu_prior_mu = 0,
    mu_prior_sigma = 100,
    dilution_stepsize = dilutionStepsize(map)
  )["mean", what]
  
  2^logtiter * 10

} 

get_titer(subset(data, sr_group == "WT-984" & ag_name == "Alpha")$titer, what = "estimate")
get_titer(subset(data, sr_group == "WT-984" & ag_name == "Alpha")$titer, what = "lower")
get_titer(subset(data, sr_group == "WT-984" & ag_name == "Alpha")$titer, what = "upper")

get_titer(subset(data, sr_group == "WT-984" & ag_name == "D614G")$titer, what = "estimate")
get_titer(subset(data, sr_group == "WT-984" & ag_name == "D614G")$titer, what = "lower")
get_titer(subset(data, sr_group == "WT-984" & ag_name == "D614G")$titer, what = "upper")

get_titer(subset(data, sr_group == "B.1+E484K" & ag_name == "Mu")$titer)
get_titer(subset(data, sr_group == "B.1+E484K" & ag_name == "B.1+E484K")$titer)

get_titer(subset(data, sr_group == "BA.5" & ag_name == "BA.5")$titer, what = "estimate")
get_titer(subset(data, sr_group == "BA.5" & ag_name == "BF.7")$titer, what = "estimate")
