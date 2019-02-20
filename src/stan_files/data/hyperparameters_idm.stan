// hyperparameter values are set to 0 if there is no prior
vector[e_K01]                  e_prior_mean01;
real                           e_prior_mean_for_intercept01;
vector[basehaz_nvars01]          e_prior_mean_for_aux01;
vector<lower=0>[e_K01]           e_prior_scale01;
real<lower=0>                  e_prior_scale_for_intercept01;
vector<lower=0>[basehaz_nvars01] e_prior_scale_for_aux01;
vector<lower=0>[e_K01]           e_prior_df01;
real<lower=0>                  e_prior_df_for_intercept01;
vector<lower=0>[basehaz_nvars01] e_prior_df_for_aux01;
real<lower=0>                  e_global_prior_scale01; // for hs priors only
real<lower=0>                  e_global_prior_df01;
real<lower=0>                  e_slab_df01;
real<lower=0>                  e_slab_scale01;

  // hyperparameter values are set to 0 if there is no prior
  vector[e_K02]                  e_prior_mean02;
  real                           e_prior_mean_for_intercept02;
  vector[basehaz_nvars02]          e_prior_mean_for_aux02;
  vector<lower=0>[e_K02]           e_prior_scale02;
  real<lower=0>                  e_prior_scale_for_intercept02;
  vector<lower=0>[basehaz_nvars02] e_prior_scale_for_aux02;
  vector<lower=0>[e_K02]           e_prior_df02;
  real<lower=0>                  e_prior_df_for_intercept02;
  vector<lower=0>[basehaz_nvars02] e_prior_df_for_aux02;
  real<lower=0>                  e_global_prior_scale02; // for hs priors only
  real<lower=0>                  e_global_prior_df02;
  real<lower=0>                  e_slab_df02;
  real<lower=0>                  e_slab_scale02;

  // hyperparameter values are set to 0 if there is no prior
  vector[e_K12]                  e_prior_mean12;
  real                           e_prior_mean_for_intercept12;
  vector[basehaz_nvars12]          e_prior_mean_for_aux12;
  vector<lower=0>[e_K12]           e_prior_scale12;
  real<lower=0>                  e_prior_scale_for_intercept12;
  vector<lower=0>[basehaz_nvars12] e_prior_scale_for_aux12;
  vector<lower=0>[e_K12]           e_prior_df12;
  real<lower=0>                  e_prior_df_for_intercept12;
  vector<lower=0>[basehaz_nvars12] e_prior_df_for_aux12;
  real<lower=0>                  e_global_prior_scale12; // for hs priors only
  real<lower=0>                  e_global_prior_df12;
  real<lower=0>                  e_slab_df12;
  real<lower=0>                  e_slab_scale12;
