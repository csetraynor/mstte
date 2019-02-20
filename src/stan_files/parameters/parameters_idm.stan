 // declares:
  //   e_{gamma,z_beta,aux_unscaled,global,local,mix,ool}
  // primitive log hazard ratios
  vector[e_K01] e_z_beta01;

  // intercept
  real e_gamma01[e_has_intercept01 == 1];

  // unscaled basehaz parameters
  //   exp model:      nvars = 0, ie. no aux parameter
  //   weibull model:  nvars = 1, ie. shape parameter
  //   gompertz model: nvars = 1, ie. scale parameter
  //   M-spline model: nvars = number of basis terms, ie. spline coefs
  //   B-spline model: nvars = number of basis terms, ie. spline coefs
  vector<lower=coefs_lb01(basehaz_type01)>[basehaz_nvars01] e_aux_unscaled01;

  // parameters for priors on log haz ratios
  real<lower=0> e_global01[e_hs01];
  vector<lower=0>[(e_hs01>0)*e_K01] e_local[e_hs01];
  real<lower=0> e_caux01[e_hs01 > 0];
  vector<lower=0>[e_K01] e_mix01[e_prior_dist01 == 5 || e_prior_dist01 == 6];
  real<lower=0> e_ool01[e_prior_dist01 == 6];
  // primitive log hazard ratios
  vector[e_K02] e_z_beta02;

  // intercept
  real e_gamma02[e_has_intercept02 == 1];

  // unscaled basehaz parameters
  //   exp model:      nvars = 0, ie. no aux parameter
  //   weibull model:  nvars = 1, ie. shape parameter
  //   gompertz model: nvars = 1, ie. scale parameter
  //   M-spline model: nvars = number of basis terms, ie. spline coefs
  //   B-spline model: nvars = number of basis terms, ie. spline coefs
  vector<lower=coefs_lb02(basehaz_type02)>[basehaz_nvars02] e_aux_unscaled02;

  // parameters for priors on log haz ratios
  real<lower=0> e_global02[e_hs02];
  vector<lower=0>[(e_hs02>0)*e_K02] e_local[e_hs02];
  real<lower=0> e_caux02[e_hs02 > 0];
  vector<lower=0>[e_K02] e_mix02[e_prior_dist02 == 5 || e_prior_dist02 == 6];
  real<lower=0> e_ool02[e_prior_dist02 == 6];

    // primitive log hazard ratios
  vector[e_K12] e_z_beta12;

  // intercept
  real e_gamma12[e_has_intercept12 == 1];

  // unscaled basehaz parameters
  //   exp model:      nvars = 0, ie. no aux parameter
  //   weibull model:  nvars = 1, ie. shape parameter
  //   gompertz model: nvars = 1, ie. scale parameter
  //   M-spline model: nvars = number of basis terms, ie. spline coefs
  //   B-spline model: nvars = number of basis terms, ie. spline coefs
  vector<lower=coefs_lb12(basehaz_type12)>[basehaz_nvars12] e_aux_unscaled12;

  // parameters for priors on log haz ratios
  real<lower=0> e_global12[e_hs12];
  vector<lower=0>[(e_hs12>0)*e_K12] e_local[e_hs12];
  real<lower=0> e_caux12[e_hs12 > 0];
  vector<lower=0>[e_K12] e_mix12[e_prior_dist12 == 5 || e_prior_dist12 == 6];
  real<lower=0> e_ool12[e_prior_dist12 == 6];
