 // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = hs
  //   4 = hs_plus
  //   5 = laplace
  //   6 = lasso
  int<lower=0,upper=6> e_prior_dist01;
  int<lower=0,upper=2> e_prior_dist_for_intercept01;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = exponential
  int<lower=0,upper=3> e_prior_dist_for_aux01; // prior for basehaz params

  // baseline hazard type:
  //   1 = weibull
  //   2 = B-splines
  //   3 = piecewise
  int<lower=1,upper=3> basehaz_type01;

  // dimensions
  int<lower=0> e_K01;                // num. predictors in event submodel
  int<lower=0> basehaz_nvars01;      // num. aux parameters for baseline hazard
  int<lower=0> qnodes01;             // num. nodes for GK quadrature
  int<lower=0> len_epts01;           // num. epts (event times)
  int<lower=0> len_qpts01;           // num. qpts (quadrature points)
  int<lower=0> len_ipts01;           // num. ipts (qpts for interval cens.)
  int<lower=0> len_cpts01;           // = len_epts + len_qpts + len_ipts
  int idx_cpts01[3,2];               // index for breaking cpts into epts,qpts,ipts

  // response and time variables
  vector[len_epts01] epts01;           // time of events
  vector[len_qpts01] qpts01;           // time at quadpoints
  vector[len_ipts01] ipts01;           // time at quadpoints for interval censoring

  // predictor matrices
  matrix[len_cpts01, e_K01] e_x01;             // predictor matrix
  vector[e_K01] e_xbar01;                    // predictor means

  // spline basis for baseline hazard
  matrix[len_epts01, basehaz_nvars01] basis_epts01;
  matrix[len_qpts01, basehaz_nvars01] basis_qpts01;
  matrix[len_ipts01, basehaz_nvars01] basis_ipts01;

  // GK quadrature weights, with (b-a)/2 scaling already incorporated
  vector[len_qpts01] qwts01;
  vector[len_ipts01] iwts01;

  // constant shift for log baseline hazard
  real norm_const01;

  // flags
  int<lower=0,upper=1> e_has_intercept01; // basehaz requires intercept

