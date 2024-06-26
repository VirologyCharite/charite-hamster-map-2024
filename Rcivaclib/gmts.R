
#' @export
srGroupGMTs <- function(map, ci_method = "HDI") {

  # Fetch titer table
  titer_table <- adjustedTiterTable(map)

  # Calculate gmts for each group and antigen
  sr_group_gmts <- lapply(
    unique(srGroups(map)),
    \(sr_group) {
      apply(
        titer_table[ , srGroups(map) == sr_group], 1, \(titers) {
          titertools::gmt(
            titers = titers,
            ci_method = ci_method,
            sigma_prior_alpha = 2,
            sigma_prior_beta = 0.75,
            mu_prior_mu = 0,
            mu_prior_sigma = 100,
            dilution_stepsize = dilutionStepsize(map)
          )['mean', 'estimate']
        }
      )
    }
  )

  # Convert to a matrix and return the values
  sr_group_gmts <- do.call(cbind, sr_group_gmts)
  colnames(sr_group_gmts) <- as.character(unique(srGroups(map)))
  sr_group_gmts <- sr_group_gmts[ ,match(levels(srGroups(map)), colnames(sr_group_gmts))]
  sr_group_gmts

}
