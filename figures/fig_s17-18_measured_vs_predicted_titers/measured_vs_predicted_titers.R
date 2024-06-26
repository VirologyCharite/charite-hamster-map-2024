
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)

source("Rcivaclib/scales.R")

set.seed(100)

map <- read.acmap("data/240123-hamster/map-continuous-fixbottom-90-corrected.ace")

# Set serum group order
srGroups(map) <- factor(
  srGroups(map),
  c("D614G sera", "Alpha sera", "Delta sera", "Beta sera",
    "B.1+E484K sera", "BA.2 sera", "BA.2-12 sera", "BA.1 sera",
    "BA.5 sera", "XBB.2 sera"
  )
)

agFillScale <- function(map) {
  fills <- agFill(map)
  names(fills) <- agNames(map)
  fills
}

#" Cross validation performed excluding 10% of titers, then predicting them, 1000 cross-validation replicates
#" were performed and 1000 optimizations were performed to find the lowest stress map each time.

crossvalidateMap <- function(
  map,
  test_proportion = 0.1,
  number_of_optimizations = 1000,
  number_of_replicates = 1000,
  optimization_number = 1,
  options = list()
) {
  
  # Perform the CV testing
  cv_results <- Racmacs:::runDimensionTestMap(
    map = map,
    dimensions_to_test = mapDimensions(map, optimization_number),
    test_proportion = test_proportion,
    minimum_column_basis = minColBasis(map, optimization_number),
    fixed_column_bases = fixedColBases(map, optimization_number),
    number_of_optimizations = number_of_optimizations,
    replicates_per_dimension = number_of_replicates,
    options = options
  )
  
  # Summarise the results
  do.call(
    bind_rows,
    lapply(seq_along(cv_results$results), \(n) {
      tibble(
        measured_titer = cv_results$titers[cv_results$results[[n]]$test_indices],
        predicted_logtiter = cv_results$results[[n]]$predictions[[1]],
        titer_index = as.vector(cv_results$results[[n]]$test_indices),
        run = n
      )
    })
  ) %>% 
    mutate(
      ag_num = Racmacs:::agNumMatrix(map)[titer_index],
      sr_num = Racmacs:::srNumMatrix(map)[titer_index]
    )
  
}

results <- crossvalidateMap(map, number_of_optimizations = 1000, number_of_replicates = 500, test_proportion = 0.1)
results %>%
  mutate(
    ag_name = agNames(map)[ag_num],
    sr_group = srGroups(map)[sr_num],
    measured_logtiter = Racmacs:::log_titers(measured_titer, 0),
    predicted_titer = as.character(round(2^predicted_logtiter*10, 6)),
    measured_titer_type = Racmacs:::titer_types_int(measured_titer),
    residual = measured_logtiter - predicted_logtiter
  ) -> results

saveRDS(results, "figures/fig_s17-18_measured_vs_predicted_titers/cv_testing_result_240123.rds")

stop()

results <- readRDS("figures/fig_s17-18_measured_vs_predicted_titers/cv_testing_result_240123.rds")

# Set detectable results subset
detectable_results <- filter(results, measured_titer_type == 1 & is.finite(residual))

# Plot histogram
detectable_results %>% 
  ggplot() + 
  geom_histogram(
    aes(
      x = residual
    ),
    binwidth = 0.2,
    fill = "darkblue",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = mean(detectable_results$residual),
    lty = 2,
    lwd = 0.6
  ) +
  coord_cartesian(
    xlim = c(-8, 8)
  ) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  labs(
    x = "Measured log titer - Predicted log titer",
    y = "Count",
    subtitle = sprintf(
      "µ = %s, σ = %s, RMSD = %s",
      round(mean(detectable_results$residual), 2),
      round(sd(detectable_results$residual), 2),
      round(sqrt(mean(detectable_results$residual^2)), 2)
    )
  ) 
ggsave("figures/fig_s17-18_measured_vs_predicted_titers/measured_vs_predicted_titers_hist_240123.png", width = 5, height = 4)

sqrt(sum(detectable_results$residual^2) / length(detectable_results$residual))

# Result without BA.1 variant and sera
detectable_results_ba1 <- subset(detectable_results, ag_name != "BA.1" & sr_group != "BA.1 convalescent")
mean(detectable_results_ba1$residual)
sd(detectable_results_ba1$residual)
sqrt(sum(detectable_results_ba1$residual^2) / length(detectable_results_ba1$residual))


# Plot split by antigen
detectable_results %>%
  mutate(
    ag_name = factor(ag_name, levels = c("D614G", "Alpha", "Delta", "Beta", "Mu", "B.1+E484K", "BA.2-12", "BA.2", "BA.1", "BA.4", "BA.5", "BF.7", "BQ.1.18", "BN.1.3.1", "XBB.2", "EG.5.1", "JN.1"))
  ) %>%
  ggplot(
    aes(
      x = ag_name,
      y = residual,
      fill = ag_name
    )
  ) +
  geom_boxplot(
    lwd = 0.25,
    outlier.shape = NA,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  scale_fill_manual(
    values = agFillScale(map)
  ) +
  titerplot_theme() +
  labs(
    x = "",
    y = "Measured log titer - Predicted log titer"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(vars(sr_group))
ggsave("figures/fig_s17-18_measured_vs_predicted_titers/measured_vs_predicted_titers_by_ag_sr_240123.png", width = 8, height=6)

