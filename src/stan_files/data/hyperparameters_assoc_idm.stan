
  // hyperparameter values are set to 0 if there is no prior
  vector[a_K01]          a_prior_mean01;
  vector<lower=0>[a_K01] a_prior_scale01;
  vector<lower=0>[a_K01] a_prior_df01;
  real<lower=0>        a_global_prior_scale01; // for hs priors only
  real<lower=0>        a_global_prior_df01;
  real<lower=0>        a_slab_df01;
  real<lower=0>        a_slab_scale01;

      // hyperparameter values are set to 0 if there is no prior
  vector[a_K02]          a_prior_mean02;
  vector<lower=0>[a_K02] a_prior_scale02;
  vector<lower=0>[a_K02] a_prior_df02;
  real<lower=0>        a_global_prior_scale02; // for hs priors only
  real<lower=0>        a_global_prior_df02;
  real<lower=0>        a_slab_df02;
  real<lower=0>        a_slab_scale02;

      // hyperparameter values are set to 0 if there is no prior
  vector[a_K12]          a_prior_mean12;
  vector<lower=0>[a_K12] a_prior_scale12;
  vector<lower=0>[a_K12] a_prior_df12;
  real<lower=0>        a_global_prior_scale12; // for hs priors only
  real<lower=0>        a_global_prior_df12;
  real<lower=0>        a_slab_df12;
  real<lower=0>        a_slab_scale12;
