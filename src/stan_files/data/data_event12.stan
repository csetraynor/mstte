 // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = hs
  //   4 = hs_plus
  //   5 = laplace
  //   6 = lasso
  int<lower=0,upper=6> e_prior_dist12;
  int<lower=0,upper=2> e_prior_dist_for_intercept12;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = exponential
  int<lower=0,upper=3> e_prior_dist_for_aux12; // prior for basehaz params

  // baseline hazard type:
  //   1 = weibull
  //   2 = B-splines
  //   3 = piecewise
  int<lower=1,upper=3> basehaz_type12;

  // dimensions
  int<lower=0> e_K12;                // num. predictors in event submodel
  int<lower=0> basehaz_nvars12;      // num. aux parameters for baseline hazard
  int<lower=0> qnodes12;             // num. nodes for GK quadrature
  int<lower=0> len_epts12;           // num. epts (event times)
  int<lower=0> len_qpts12;           // num. qpts (quadrature points)
  int<lower=0> len_ipts12;           // num. ipts (qpts for interval cens.)
  int<lower=0> len_cpts12;           // = len_epts + len_qpts + len_ipts
  int idx_cpts12[3,2];               // index for breaking cpts into epts,qpts,ipts

  // response and time variables
  vector[len_epts12] epts12;           // time of events
  vector[len_qpts12] qpts12;           // time at quadpoints
  vector[len_ipts12] ipts12;           // time at quadpoints for interval censoring

  // predictor matrices
  matrix[len_cpts12, e_K12] e_x12;             // predictor matrix
  vector[e_K12] e_xbar12;                    // predictor means

  // spline basis for baseline hazard
  matrix[len_epts12, basehaz_nvars12] basis_epts12;
  matrix[len_qpts12, basehaz_nvars12] basis_qpts12;
  matrix[len_ipts12, basehaz_nvars12] basis_ipts12;

  // GK quadrature weights, with (b-a)/2 scaling already incorporated
  vector[len_qpts12] qwts12;
  vector[len_ipts12] iwts12;

  // constant shift for log baseline hazard
  real norm_const12;

  // flags
  int<lower=0,upper=1> e_has_intercept12; // basehaz requires intercept

