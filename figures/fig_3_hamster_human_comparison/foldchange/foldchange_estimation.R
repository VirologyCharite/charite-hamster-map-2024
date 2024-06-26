
rm(list = ls())

library(Racmacs)
library(cmdstanr)
library(bayesplot)
library(tidyverse)

library(titertools)

source('Rcivaclib/gmts.R')

# Read in data
hamster = read.acmap('data/240123-hamster/map-continuous-fixbottom-90-corrected.ace')
human = read.acmap('data/220511/map-full-PRNT-220511.ace')

# Re-name human data so it's the same as hamster data
agNames(human)[agNames(human) == "B.1.351"] <- 'Beta'
agNames(human)[agNames(human) == "B.1.1.7"] <- 'Alpha'
agNames(human)[agNames(human) == "B.1.617.2"] <- 'Delta'

srGroups(human) <- recode_factor(srGroups(human), "B.1.351 convalescent" = "Beta sera")
srGroups(human) <- recode_factor(srGroups(human), "WT convalescent" = "D614G sera")
srGroups(human) <- recode_factor(srGroups(human), "B.1.1.7 convalescent" = "Alpha sera")

# Subset maps to common antigens and sera
hamster_subset = subsetMap(hamster, 
                           antigens = c("D614G", "Alpha", "Beta", "Delta"),
                           sera = srNames(hamster)[srGroups(hamster) %in% c("D614G sera", "Alpha sera", "Beta sera")])

human_subset = subsetMap(human,
                         antigens = c("D614G", "Alpha", "Beta", "Delta"),
                         sera = srNames(human)[srGroups(human) %in% c("D614G sera", "Alpha sera", "Beta sera")])

mapName(hamster_subset) <- 'hamster'
mapName(human_subset) <- 'human'

# Combine data
data = mergeMaps(hamster_subset, human_subset)

# Homologous ag info
homologous_ags <- c(
  "Beta sera hamster" = "Beta",
  "D614G sera hamster" = "D614G",
  "Alpha sera hamster" = "Alpha",
  "Beta sera human" = "Beta",
  "D614G sera human" = "D614G",
  "Alpha sera human" = "Alpha"
)

slope_sr_groups <- c('D614G sera', 'Alpha sera', 'Beta sera')

# Arrange the data
srGroups(data) <- factor(srGroups(data), levels = slope_sr_groups)

# Run the model
options(mc.cores = parallel::detectCores())

# Build up titer information
titer_layers <- Racmacs:::titerTableLayers(data)

upper_logtiter_lims <- Racmacs:::logtiterTableLayers(data)
lower_logtiter_lims <- Racmacs:::logtiterTableLayers(data)

for (map_ in layerNames(data)) {
  for (serum_group in slope_sr_groups) {
    # get homologous ag
    homologous_ag <- unname(homologous_ags[paste(serum_group, map_)])
    
    for (ag in agNames(data)) {
      diff_lims <- titertools:::calc_titer_diff_lims(
        titers1 = titer_layers[[map_]][match(homologous_ag, agNames(data)), which(srGroups(data) %in% c(serum_group))],
        titers2 = titer_layers[[map_]][match(ag, agNames(data)), which(srGroups(data) %in% c(serum_group))],
        dilution_stepsize = 0
      )
      upper_logtiter_lims[[map_]][match(ag, agNames(data)), which(srGroups(data) %in% c(serum_group))] <- diff_lims$max_diffs
      lower_logtiter_lims[[map_]][match(ag, agNames(data)), which(srGroups(data) %in% c(serum_group))] <- diff_lims$min_diffs
    }
  }
}

# Replace NAs with 0s
upper_logtiter_lims <- lapply(upper_logtiter_lims, \(upper_lims) {
  upper_lims[is.na(upper_lims)] <- 0
  upper_lims
})

lower_logtiter_lims <- lapply(lower_logtiter_lims, \(lower_lims) {
  lower_lims[is.na(lower_lims)] <- 0
  lower_lims
})


# Assemble the data list for the model
data_list <- list(
  only_prior = 0,
  N_datasets = numLayers(data),
  N_ags = numAntigens(data),
  N_srs = numSera(data),
  N_sr_groups = numSeraGroups(data),
  sr_groups = as.numeric(srGroups(data)),
  titertypes = simplify2array(Racmacs:::titertypesTableLayers(data)),
  lower_logtiter_lims = simplify2array(lower_logtiter_lims),
  upper_logtiter_lims = simplify2array(upper_logtiter_lims)
)

# Decide on starting conditions for the model
data_init <- list(
  logtiter_error_sigma = rep(0.8, numLayers(data)),
  dataset_slope_effect = rep(1, numLayers(data))
)

# Fetch the model
mod <- cmdstan_model(stan_file = 'figures/hamster_human_comparison/foldchange/slope_calculation.stan')

# Sample
mod_optimized <- mod$sample(
  data = data_list,
  init = list(data_init, data_init, data_init, data_init),
  chains = 4, 
  refresh = 500,
  iter_sampling  = 10000,
  iter_warmup = 1000
)

# Draw from the posterior, map slope effect
mod_optimized_draws <- mod_optimized$draws(format = "df", variables = c('dataset_slope_effect'))

# Re-name columns
colnames(mod_optimized_draws) <- c("hamster", "human", ".chain", ".iteration", ".draw")

# Save the draws
saveRDS(mod_optimized_draws, 'figures/hamster_human_comparison/foldchange/slope_calculation_dataset_slope_effect_draws.rds')

# Draw from the posterior, ag_folddrops
mod_optimized_draws_foldchange <- mod_optimized$draws(format = "df", variables = c('ag_folddrops'))

# Save the draws
saveRDS(mod_optimized_draws_foldchange, 'figures/hamster_human_comparison/foldchange/slope_calculation_ag_folddrops_draws.rds')

# Check convergence
summary <- mod_optimized$summary()
summary
# # A tibble: 17 × 10
# variable                      mean     median     sd    mad       q5      q95  rhat ess_bulk ess_tail
# <chr>                        <dbl>      <dbl>  <dbl>  <dbl>    <dbl>    <dbl> <dbl>    <dbl>    <dbl>
# 1 lp__                    -400.      -400.      3.77   3.70   -407.    -394.     1.00    5750.   12493.
# 2 dataset_slope_effect[1]    0.569      0.544   0.145  0.125     0.386    0.840  1.00    4856.    7891.
# 3 dataset_slope_effect[2]    0.823      0.787   0.199  0.169     0.572    1.19   1.00    4490.    7392.
# 4 ag_folddrops[1,1]         -0.00356   -0.00328 0.162  0.152    -0.269    0.259  1.00   41574.   24323.
# 5 ag_folddrops[2,1]          0.716      0.696   0.221  0.217     0.386    1.11   1.00    7907.   16376.
# 6 ag_folddrops[3,1]         -2.59      -2.56    0.589  0.590    -3.59    -1.67   1.00    4817.    7974.
# 7 ag_folddrops[4,1]         -0.605     -0.581   0.782  0.725    -1.92     0.625  1.00   29345.   20546.
# 8 ag_folddrops[1,2]         -3.77      -3.73    0.867  0.865    -5.26    -2.41   1.00    5113.    8369.
# 9 ag_folddrops[2,2]         -0.0121    -0.0116  0.338  0.319    -0.569    0.540  1.00   40595.   23982.
# 10 ag_folddrops[3,2]         -4.43      -4.39    0.996  1.00     -6.14    -2.86   1.00    4901.    8281.
# 11 ag_folddrops[4,2]         -5.40      -5.37    1.20   1.20     -7.45    -3.51   1.00    4733.    7989.
# 12 ag_folddrops[1,3]         -4.47      -4.43    0.989  0.987    -6.16    -2.91   1.00    4670.    7780.
# 13 ag_folddrops[2,3]         -0.461     -0.442   0.216  0.211    -0.844   -0.138  1.00   14667.   19311.
# 14 ag_folddrops[3,3]         -0.00498   -0.00318 0.194  0.182    -0.324    0.312  1.00   37058.   22953.
# 15 ag_folddrops[4,3]         -4.30      -4.26    0.958  0.961    -5.93    -2.79   1.00    4729.    7880.
# 16 logtiter_error_sigma[1]    0.723      0.713   0.107  0.102     0.569    0.914  1.00   33670.   25576.
# 17 logtiter_error_sigma[2]    0.889      0.888   0.0417 0.0412    0.824    0.961  1.00   39110.   29182.

# Look at traces
gp <- mcmc_trace(mod_optimized$draws(variables = c('dataset_slope_effect')))
gp
ggsave('figures/hamster_human_comparison/foldchange/slope_calculation_dataset_slope_effect_traceplots.png', gp, width = 15, height = 10)
