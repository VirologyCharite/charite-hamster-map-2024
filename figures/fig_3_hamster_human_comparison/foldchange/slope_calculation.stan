

functions {
  
  // Function to calculate likelihood of interval censored measurement given 
  // mean and standard deviation of a normal distribution
  real normal_int_censored_likelihood(
    real lower_lim, 
    real upper_lim, 
    real mu, 
    real sigma
    ) {
      
      real result;
      
      if (is_inf(lower_lim) && is_inf(upper_lim)) {
        
        result = 0;
        
      } else if (lower_lim == upper_lim) {
        
        result = normal_lpdf(
          lower_lim | mu, sigma
          );
          
      } else if (!is_inf(lower_lim) && !is_inf(upper_lim)) {
        
        result = log_diff_exp(
          normal_lcdf(upper_lim | mu, sigma),
          normal_lcdf(lower_lim | mu, sigma)  // shouldn't this be normal_lccdf?
          );
          
      } else if (!is_inf(lower_lim) && is_inf(upper_lim)) {
        
        result = normal_lccdf(
          lower_lim | mu, sigma
          );
          
      } else if (is_inf(lower_lim) && !is_inf(upper_lim)) {
        
        result = normal_lcdf(
          upper_lim | mu, sigma
          );
          
      }
      
      return result;
      
    }
    
}


data {
  
  // Flag for sampling from the prior alone
  int only_prior;
  
  // Dataset related data
  int<lower=0> N_datasets;
  
  // Antigen related data
  int<lower=0> N_ags;
  
  // Sera related data
  int<lower=0> N_srs;
  int<lower=0> N_sr_groups;
  int sr_groups[N_srs];
  
  // Titer-measurement related data
  real titertypes[N_ags, N_srs, N_datasets];
  real lower_logtiter_lims[N_ags, N_srs, N_datasets];
  real upper_logtiter_lims[N_ags, N_srs, N_datasets];
  
}

parameters {
  
  // Slope related parameters
  vector<lower=-0.05>[N_datasets] dataset_slope_effect;
  
  // Antigen folddrop related parameters
  matrix[N_ags, N_sr_groups] ag_folddrops;
  
  // Variation related data
  vector<lower=0>[N_datasets] logtiter_error_sigma;
  
}

model {
  
  // Prior for the dataset slope effect
  dataset_slope_effect ~ normal(1, 10);
  
  // Prior for the log titer error sigma
  logtiter_error_sigma ~ inv_gamma(2, 5);
  
  // Prior for ag_folddrops
  for (ag in 1:N_ags) {
    ag_folddrops[ag] ~ normal(-1, 3);
  }
  
  // Cycle through each antigen and sera to get the predicted foldchange
  // and likelihood of the measurements given that prediction
  if (only_prior == 0) {
    for (dataset in 1:N_datasets) {
      for (sr in 1:N_srs) {
        for (ag in 1:N_ags) {
          
          if (titertypes[ag, sr, dataset] > 0) {
            
            // Fetch serum group number
            int sr_group = sr_groups[sr];
            
            // Predict foldchange
            real predicted_foldchange = ag_folddrops[ag, sr_group] * dataset_slope_effect[dataset];
            
            // Accumulate likelihood of predicted foldchange
            target += normal_int_censored_likelihood(
              lower_logtiter_lims[ag, sr, dataset],
              upper_logtiter_lims[ag, sr, dataset],
              predicted_foldchange,
              logtiter_error_sigma[dataset]
              );
          }
        }
      }
    }
  }
}

