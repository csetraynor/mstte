 // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = hs
  //   4 = hs_plus
  //   5 = laplace
  //   6 = lasso
  int<lower=0,upper=6> e_prior_dist02;
  int<lower=0,upper=2> e_prior_dist_for_intercept02;

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = exponential
  int<lower=0,upper=3> e_prior_dist_for_aux02; // prior for basehaz params

  // baseline hazard type:
  //   1 = weibull
  //   2 = B-splines
  //   3 = piecewise
  int<lower=1,upper=3> basehaz_type02;

  // dimensions
  int<lower=0> e_K02;                // num. predictors in event submodel
  int<lower=0> basehaz_nvars02;      // num. aux parameters for baseline hazard
  int<lower=0> qnodes02;             // num. nodes for GK quadrature
  int<lower=0> len_epts02;           // num. epts (event times)
  int<lower=0> len_qpts02;           // num. qpts (quadrature points)
  int<lower=0> len_ipts02;           // num. ipts (qpts for interval cens.)
  int<lower=0> len_cpts02;           // = len_epts + len_qpts + len_ipts
  int idx_cpts02[3,2];               // index for breaking cpts into epts,qpts,ipts

  // response and time variables
  vector[len_epts02] epts02;           // time of events
  vector[len_qpts02] qpts02;           // time at quadpoints
  vector[len_ipts02] ipts02;           // time at quadpoints for interval censoring

  // predictor matrices
  matrix[len_cpts02, e_K02] e_x02;             // predictor matrix
  vector[e_K02] e_xbar02;                    // predictor means

  // spline basis for baseline hazard
  matrix[len_epts02, basehaz_nvars02] basis_epts02;
  matrix[len_qpts02, basehaz_nvars02] basis_qpts02;
  matrix[len_ipts02, basehaz_nvars02] basis_ipts02;

  // GK quadrature weights, with (b-a)/2 scaling already incorporated
  vector[len_qpts02] qwts02;
  vector[len_ipts02] iwts02;

  // constant shift for log baseline hazard
  real norm_const02;

  // flags
  int<lower=0,upper=1> e_has_intercept02; // basehaz requires intercept

