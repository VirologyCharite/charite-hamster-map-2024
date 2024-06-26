

calculate_fold_change <- function(map_info, sr_groups, homologous_info) {
  foldchange_table <- tibble(
    ag_name = character(0),
    sr_group = character(0),
    map = character(0),
    mean_diff = numeric(0),
    mean_diff_lower = numeric(0),
    mean_diff_upper = numeric(0)
  )
  
  # Compute fold drops for map
  for (sg in sr_groups) {
    for (map_name in unique(map_info$map)) {
      homologous_ag <- unname(homologous_ags[paste(sg, map_name)])
      homologous_titers <- subset(map_info, sr_group == sg & ag_name == homologous_ag & map == map_name)$titer
        
      for (ag in unique(subset(map_info, sr_group == sg & map == map_name)$ag_name)) {
        print(paste(sg, map_name, ag, homologous_ag))
        if (ag == homologous_ag) {
          if (length(setdiff(unique(homologous_titers), c("."))) == 0) {
            # The homologous antigen didn't get titrated, add NA
            print('is homologous and doesnt exist difference NA')
            foldchange_table <- bind_rows(
              foldchange_table,
              tibble(
                ag_name = ag,
                sr_group = sg,
                map = map_name,
                mean_diff = NA,
                mean_diff_lower = NA,
                mean_diff_upper = NA
              )
            )
          } else {
            # The homologous ag did get titrated, the difference will be 0, therefore, add mean_diff = 0.
            print('is homologous, exists, difference 0')
            foldchange_table <- bind_rows(
              foldchange_table,
              tibble(
                ag_name = ag,
                sr_group = sg,
                map = map_name,
                mean_diff = 0,
                mean_diff_lower = NA,
                mean_diff_upper = NA
              )
            )
          }
        }
        if (ag != homologous_ag) {
          titers1 <- subset(map_info, sr_group == sg & ag_name == ag & map == map_name)$titer

          titers <- tibble(titers1, homologous_titers)
          titers <- subset(titers, !(titers1 %in% c(".", "*")) & !(homologous_titers %in% c(".", "*")))
          print(titers$homologous_titers)
          print(titers$titers1)

          if (length(titers$homologous_titers) >= 1) {
            # Both the homologous and the comparison antigen got titrated
            print('Both the homologous and the comparison antigen got titrated, calculating')
            fold_drop <- titertools::log2diff(
              titers2 = titers$titers1,
              titers1 = titers$homologous_titers,
              dilution_stepsize = 0,
              ci_method = "HDI",
              mu_prior_mu = 0,
              mu_prior_sigma = 100,
              sigma_prior_alpha = 2,
              sigma_prior_beta = 0.75
            )

            foldchange_table <- bind_rows(
              foldchange_table,
              tibble(
                ag_name = ag,
                sr_group = sg,
                map = map_name,
                mean_diff = as.numeric(fold_drop[1]),
                mean_diff_lower = as.numeric(fold_drop[3]),
                mean_diff_upper = as.numeric(fold_drop[5])
              )
            )
          } else {
            # The comparison antigen did not get titrated
            print('The comparison antigen did not get titrated, adding NA')
            foldchange_table <- bind_rows(
              foldchange_table,
              tibble(
                ag_name = ag,
                sr_group = sg,
                map = map_name,
                mean_diff = NA,
                mean_diff_lower = NA,
                mean_diff_upper = NA
              )
            )
          }
        }
      }
    }
  }
  foldchange_table
}


