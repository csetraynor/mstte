  // declares:
  //   a_{z_beta,global,local,mix,ool}
  vector[a_K01] a_z_beta01; // primitive assoc params

  // parameters for priors on assoc params
  real<lower=0> a_global01[a_hs01];
  vector<lower=0>[(a_hs01>0)*a_K01] a_local01[a_hs01];
  real<lower=0> a_caux01[a_hs01 > 0];
  vector<lower=0>[a_K01] a_mix01[a_prior_dist01 == 5 || a_prior_dist01 == 6];
  real<lower=0> a_ool01[a_prior_dist01 == 6];

    //   a_{z_beta,global,local,mix,ool}
  vector[a_K02] a_z_beta02; // primitive assoc params

  // parameters for priors on assoc params
  real<lower=0> a_global02[a_hs02];
  vector<lower=0>[(a_hs02>0)*a_K02] a_local02[a_hs02];
  real<lower=0> a_caux02[a_hs02 > 0];
  vector<lower=0>[a_K02] a_mix02[a_prior_dist02 == 5 || a_prior_dist02 == 6];
  real<lower=0> a_ool02[a_prior_dist02 == 6];

    //   a_{z_beta,global,local,mix,ool}
  vector[a_K12] a_z_beta12; // primitive assoc params

  // parameters for priors on assoc params
  real<lower=0> a_global12[a_hs12];
  vector<lower=0>[(a_hs12>0)*a_K12] a_local12[a_hs12];
  real<lower=0> a_caux12[a_hs12 > 0];
  vector<lower=0>[a_K12] a_mix12[a_prior_dist12 == 5 || a_prior_dist12 == 6];
  real<lower=0> a_ool12[a_prior_dist12 == 6];
