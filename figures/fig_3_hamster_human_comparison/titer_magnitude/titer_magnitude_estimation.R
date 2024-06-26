
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

# Combine data
data = mergeMaps(hamster_subset, human_subset)

# Run the model
options(mc.cores = parallel::detectCores())

# Replace logtiter table layer NAs with 0s, they will be dealt with by titer types
logtiter_layers <- Racmacs:::logtiterTableLayers(data)
logtiter_layers <- lapply(logtiter_layers, \(logtiters) {
  logtiters[is.na(logtiters)] <- 0
  logtiters
})

# Assemble the data list for the model
data_list <- list(
  only_prior = 0,
  N_datasets = numLayers(data),
  N_ags = numAntigens(data),
  N_srs = numSera(data),
  N_sr_groups = numSeraGroups(data),
  sr_groups = as.numeric(srGroups(data)),
  logtiters = simplify2array(logtiter_layers),
  titertypes = simplify2array(Racmacs:::titertypesTableLayers(data))
)

# Decide on starting conditions for the model
sr_group_gmts <- srGroupGMTs(data, ci_method = "quap")
sr_group_gmts[is.na(sr_group_gmts)] <- 7

data_init <- list(
  sr_group_gmts = sr_group_gmts,
  sr_individual_effects = rep(0, numSera(data)),
  map_effects = rep(0, numLayers(data)),
  logtiter_error_sigma = rep(0.8, numLayers(data))
)

# Fetch the model
mod <- cmdstan_model("figures/hamster_human_comparison/titer_magnitude/titer_magnitude.stan")

# Run the model
mod_optimized <- mod$sample(
  data = data_list,
  init = list(data_init, data_init, data_init, data_init),
  chains = 4, 
  refresh = 500,
  iter_sampling  = 10000,
  iter_warmup = 1000
)

# Draw from the posterior
draws <- mod_optimized$draws(format = "df", variables = c('dataset_effects'))

# Re-name columns
colnames(draws) <- c("hamster", "human", ".chain", ".iteration", ".draw")

# Save the draws
saveRDS(draws, 'figures/hamster_human_comparison/titer_magnitude/dataset_magnitude_effect_posterior_samples.rds')

# Check convergence
mod_optimized$summary("dataset_effects")
# # A tibble: 2 × 10
# variable            mean median    sd   mad    q5   q95  rhat ess_bulk ess_tail
# <chr>              <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
# 1 dataset_effects[1]  0.952  0.928  3.66  3.67 -5.05  7.00  1.00    2127.    4233.
# 2 dataset_effects[2] -2.88  -2.88   3.47  3.45 -8.67  2.87  1.00    1707.    3055.

# Look at traces
mcmc_trace(mod_optimized$draws(variables = c('dataset_effects'))) -> gp
gp
ggsave('figures/hamster_human_comparison/titer_magnitude/dataset_magnitude_effect_posterior_samples_traces.png', width = 15, height = 12)
