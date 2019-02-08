functions {

  /**
    * Log hazard for exponential distribution
  *
    * @param eta Vector, linear predictor
  * @return A vector
  */
    vector exponential_log_haz(vector eta) {
      return eta;
    }

  /**
    * Log hazard for Weibull distribution
  *
    * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param shape Real, Weibull shape
  * @return A vector
  */
    vector weibull_log_haz(vector eta, vector t, real shape) {
      vector[rows(eta)] res;
      res = log(shape) + (shape - 1) * log(t) + eta;
      return res;
    }

  /**
    * Log hazard for Gompertz distribution
  *
    * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param scale Real, Gompertz scale
  * @return A vector
  */
    vector gompertz_log_haz(vector eta, vector t, real scale) {
      vector[rows(eta)] res;
      res = scale * t + eta;
      return res;
    }

  /**
    * Log hazard for B-spline model
  *
    * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param coefs Vector, B-spline coefficients
  * @return A vector
  */
    vector bspline_log_haz(vector eta, matrix basis, vector coefs) {
      vector[rows(eta)] res;
      res = basis * coefs + eta;
      return res;
    }

  /**
    * Evaluate log survival or log CDF from the log hazard evaluated at
  * quadrature points and a corresponding vector of quadrature weights
  *
    * @param qwts Vector, the quadrature weights
  * @param log_hazard Vector, log hazard at the quadrature points
  * @param qnodes Integer, the number of quadrature points for each individual
  * @param N Integer, the number of individuals (ie. rows(log_hazard) / qnodes)
  * @return A vector
  */
    real quadrature_log_surv(vector qwts, vector log_hazard) {
      real res;
      res = - dot_product(qwts, exp(log_hazard)); // sum across all individuals
      return res;
    }

  vector quadrature_log_cdf(vector qwts, vector log_hazard, int qnodes, int N) {
    int M = rows(log_hazard);
    vector[M] hazard = exp(log_hazard);
    matrix[N,qnodes] qwts_mat = to_matrix(qwts,   N, qnodes);
    matrix[N,qnodes] haz_mat  = to_matrix(hazard, N, qnodes);
    vector[N] chaz = rows_dot_product(qwts_mat, haz_mat);
    vector[N] res;
    res = log(1 - exp(- chaz));
    return res;
  }

  vector quadrature_log_cdf2(vector qwts_lower, vector log_hazard_lower,
                             vector qwts_upper, vector log_hazard_upper,
                             int qnodes, int N) {
    int M = rows(log_hazard_lower);
    vector[M] hazard_lower = exp(log_hazard_lower);
    vector[M] hazard_upper = exp(log_hazard_upper);
    matrix[N,qnodes] qwts_lower_mat = to_matrix(qwts_lower,   N, qnodes);
    matrix[N,qnodes] qwts_upper_mat = to_matrix(qwts_upper,   N, qnodes);
    matrix[N,qnodes] haz_lower_mat  = to_matrix(hazard_lower, N, qnodes);
    matrix[N,qnodes] haz_upper_mat  = to_matrix(hazard_upper, N, qnodes);
    vector[N] chaz_lower = rows_dot_product(qwts_lower_mat, haz_lower_mat);
    vector[N] chaz_upper = rows_dot_product(qwts_upper_mat, haz_upper_mat);
    vector[N] surv_lower = exp(- chaz_lower);
    vector[N] surv_upper = exp(- chaz_upper);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }


  /**
    * Log hazard for M-spline model
  *
    * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param coefs Vector, M-spline coefficients
  * @return A vector
  */
    vector mspline_log_haz(vector eta, matrix basis, vector coefs) {
      vector[rows(eta)] res;
      res = log(basis * coefs) + eta;
      return res;
    }

  /**
    * Cornish-Fisher expansion for standard normal to Student t
  *
    * See result 26.7.5 of
  * http://people.math.sfu.ca/~cbm/aands/page_949.htm
  *
    * @param z A scalar distributed standard normal
  * @param df A scalar degrees of freedom
  * @return An (approximate) Student t variate with df degrees of freedom
  */
    real CFt(real z, real df) {
      real z2 = square(z);
      real z3 = z2 * z;
      real z5 = z2 * z3;
      real z7 = z2 * z5;
      real z9 = z2 * z7;
      real df2 = square(df);
      real df3 = df2 * df;
      real df4 = df2 * df2;
      return z + (z3 + z) / (4 * df) + (5 * z5 + 16 * z3 + 3 * z) / (96 * df2)
      + (3 * z7 + 19 * z5 + 17 * z3 - 15 * z) / (384 * df3)
      + (79 * z9 + 776 * z7 + 1482 * z5 - 1920 * z3 - 945 * z) / (92160 * df4);
    }

  /**
    * Return the lower bound for the baseline hazard parameters
  *
    * @param type An integer indicating the type of baseline hazard
  * @return A real
  */
    real coefs_lb(int type) {
      real lb;
      if (type == 2) // B-splines, on log haz scale
      lb = negative_infinity();
      else if (type == 3) // piecewise constant, on log haz scale
      lb = negative_infinity();
      else
        lb = 0;
      return lb;
    }
 /**
  * Scale the primitive population level parameters based on prior information
  *
  * @param z_beta A vector of primitive parameters
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Vectors of mean and scale parameters
  *   for the prior distributions
  * @return A vector containing the population level parameters (coefficients)
  */
  vector make_beta(vector z_beta, int prior_dist, vector prior_mean,
                   vector prior_scale, vector prior_df) {
    vector[rows(z_beta)] beta;
    if (prior_dist == 0) beta = z_beta;
    else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
    else if (prior_dist == 2) for (k in 1:rows(prior_mean)) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
    }
    return beta;
  }


  /**
    * Log-prior for coefficients
  *
    * @param z_beta Vector of primative coefficients
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_scale Real, scale for the prior distribution
  * @param prior_df Real, df for the prior distribution
  * @param global_prior_df Real, df for the prior for the global hs parameter
  * @param local Vector of hs local parameters
  * @param global Real, the global parameter
  * @param mix Vector of shrinkage parameters
  * @param one_over_lambda Real
  * @return Nothing
  */
    real beta_lp(vector z_beta, int prior_dist, vector prior_scale,
               vector prior_df ) {
    if      (prior_dist == 1) target += normal_lpdf(z_beta | 0, 1);
    else if (prior_dist == 2) target += normal_lpdf(z_beta | 0, 1); // Student t
        /* else prior_dist is 0 and nothing is added */
          return target();
    }

  /**
  * Log-prior for baseline hazard parameters
  *
  * @param aux_unscaled Vector (potentially of length 1) of unscaled
  *   auxiliary parameter(s)
  * @param dist Integer specifying the type of prior distribution
  * @param df Real specifying the df for the prior distribution
  * @return Nothing
  */
  real basehaz_lp(vector aux_unscaled, int dist, vector df) {
    if (dist > 0) {
      if (dist == 1)
        target += normal_lpdf(aux_unscaled | 0, 1);
      else if (dist == 2)
        target += student_t_lpdf(aux_unscaled | df, 0, 1);
      else
        target += exponential_lpdf(aux_unscaled | 1);
    }
    return target();
  }

  /**
  * Log-prior for intercept parameters
  *
  * @param gamma Real, the intercept parameter
  * @param dist Integer, the type of prior distribution
  * @param mean Real, mean of prior distribution
  * @param scale Real, scale for the prior distribution
  * @param df Real, df for the prior distribution
  * @return Nothing
  */
  real gamma_lp(real gamma, int dist, real mean, real scale, real df) {
    if (dist == 1)  // normal
      target += normal_lpdf(gamma | mean, scale);
    else if (dist == 2)  // student_t
      target += student_t_lpdf(gamma | df, mean, scale);
    /* else dist is 0 and nothing is added */
    return target();
  }

  /**
    * Raise each element of x to the power of y
  *
    * @param x Vector
  * @param y Real, the power to raise to
  * @return vector
  */
    vector pow_vec(vector x, real y) {
      int N = rows(x);
      vector[N] res;
      for (n in 1:N)
        res[n] = pow(x[n], y);
      return res;
    }

  /**
    * Log survival and log CDF for exponential distribution
  *
    * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @return A vector
  */
    vector exponential_log_surv(vector eta, vector t) {
      vector[rows(eta)] res;
      res = - t .* exp(eta);
      return res;
    }

  vector exponential_log_cdf(vector eta, vector t) {
    vector[rows(eta)] res;
    res = log(1 - exp(-t .* exp(eta)));
    return res;
  }

  vector exponential_log_cdf2(vector eta, vector t_lower, vector t_upper) {
    int N = rows(eta);
    vector[N] exp_eta = exp(eta);
    vector[N] surv_lower = exp(-t_lower .* exp_eta);
    vector[N] surv_upper = exp(-t_upper .* exp_eta);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }

  /**
    * Log survival and log CDF for Weibull distribution
  *
    * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param shape Real, Weibull shape
  * @return A vector
  */
    vector weibull_log_surv(vector eta, vector t, real shape) {
      vector[rows(eta)] res;
      res = - pow_vec(t, shape) .* exp(eta);
      return res;
    }

  vector weibull_log_cdf(vector eta, vector t, real shape) {
    vector[rows(eta)] res;
    res = log(1 - exp(- pow_vec(t, shape) .* exp(eta)));
    return res;
  }

  vector weibull_log_cdf2(vector eta, vector t_lower, vector t_upper, real shape) {
    int N = rows(eta);
    vector[N] exp_eta = exp(eta);
    vector[N] surv_lower = exp(- pow_vec(t_lower, shape) .* exp_eta);
    vector[N] surv_upper = exp(- pow_vec(t_upper, shape) .* exp_eta);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }

  /**
    * Log survival and log CDF for Gompertz distribution
  *
    * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param scale Real, Gompertz scale
  * @return A vector
  */
    vector gompertz_log_surv(vector eta, vector t, real scale) {
      vector[rows(eta)] res;
      res = inv(scale) * -(exp(scale * t) - 1) .* exp(eta);
      return res;
    }

  vector gompertz_log_cdf(vector eta, vector t, real scale) {
    vector[rows(eta)] res;
    res = log(1 - exp(inv(scale) * -(exp(scale * t) - 1) .* exp(eta)));
    return res;
  }

  vector gompertz_log_cdf2(vector eta, vector t_lower, vector t_upper, real scale) {
    int N = rows(eta);
    real inv_scale = inv(scale);
    vector[N] exp_eta = exp(eta);
    vector[N] surv_lower = exp(inv_scale * -(exp(scale * t_lower) - 1) .* exp_eta);
    vector[N] surv_upper = exp(inv_scale * -(exp(scale * t_upper) - 1) .* exp_eta);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }

  /**
    * Log survival and log CDF for M-spline model
  *
    * @param eta Vector, linear predictor
  * @param t Vector, event or censoring times
  * @param coefs Vector, M-spline coefficients
  * @return A vector
  */
    vector mspline_log_surv(vector eta, matrix ibasis, vector coefs) {
      vector[rows(eta)] res;
      res = - (ibasis * coefs) .* exp(eta);
      return res;
    }

  vector mspline_log_cdf(vector eta, matrix ibasis, vector coefs) {
    vector[rows(eta)] res;
    res = log(1 - exp(-(ibasis * coefs) .* exp(eta)));
    return res;
  }

  vector mspline_log_cdf2(vector eta, matrix ibasis_lower, matrix ibasis_upper, vector coefs) {
    int N = rows(eta);
    vector[N] exp_eta = exp(eta);
    vector[N] surv_lower = exp(-(ibasis_lower * coefs) .* exp_eta);
    vector[N] surv_upper = exp(-(ibasis_upper * coefs) .* exp_eta);
    vector[N] res;
    res = log(surv_lower - surv_upper);
    return res;
  }

}

data {
  int<lower=1> nt;         // number of transitions
  int<lower=0> s_K[nt];      // size of fixed effects for each transition
  int<lower=0> nK;         // total number of fixed effects
  int<lower=0> s_vars[nt];  // basis for each transition
  int<lower=0> Nvars;      // total number of vars for basis
  int<lower=0> s_event[nt]; // events for each transition
  int<lower=0> s_rcens[nt]; // censoreds for each transition
  int<lower=0> Nevent;     // total number of events
  int<lower=0> Nrcens;     // total number of right censored

  // log crude event rate (used for centering log baseline hazard)
  vector[nt] log_crude_event_rate;

  // response and time variables
  vector[nK] x_bar;           // predictor means

  vector[Nevent] t_event;  // time of events
  vector[Nrcens] t_rcens;  // time of right censoring

  int<lower=0> Nxevent;
  int<lower=0> Nxrcens;
  vector[Nxevent] x_event;
  vector[Nxrcens] x_rcens;

  // num. aux parameters for baseline hazard
  int<lower=0> Nbasis_event;
  int<lower=0> Nibasis_event;
  int<lower=0> Nibasis_rcens;

  vector[Nbasis_event] basis_event;
  vector[Nibasis_event] ibasis_event;
  vector[Nibasis_rcens] ibasis_rcens;

  // baseline hazard type:
  //   1 = weibull
  //   2 = B-splines
  //   3 = piecewise
  //   4 = M-splines
  //   5 = exponential
  //   6 = gompertz
  int<lower=1,upper=7> type[nt];

  // flags
  int<lower=0,upper=1> has_intercept[nt]; // basehaz requires intercept
  int<lower=0> N_has_intercept;
  int<lower=0,upper=1> prior_PD;      // draw only from prior predictive dist.

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  int<lower=0,upper=2> prior_dist[nt];

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  int<lower=0,upper=2> prior_dist_for_intercept[nt];

  // prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = exponential
  int<lower=0,upper=3> prior_dist_for_aux[nt];


  // hyperparameter (log hazard ratios), set to 0 if there is no prior
  vector[nK] prior_mean;
  vector<lower=0>[nK] prior_scale;

  vector<lower=0>[nK]  prior_df;
  real<lower=0>       global_prior_scale[nt]; // for hs priors only
  real<lower=0>       global_prior_df[nt];
  real<lower=0>       slab_scale[nt];
  real<lower=0>       slab_df[nt];

  // hyperparameters (intercept), set to 0 if there is no prior
  real                prior_mean_for_intercept[nt];
  real<lower=0>       prior_scale_for_intercept[nt];
  real<lower=0>       prior_df_for_intercept[nt];

  // hyperparameters (basehaz pars), set to 0 if there is no prior
  vector<lower=0>[Nvars] prior_scale_for_aux;
  vector<lower=0>[Nvars] prior_df_for_aux;
}

parameters {
  // log hazard ratios
  vector[nK] z_beta;

  // unscaled basehaz parameters
  //   exp model:      nvars = 0, ie. no aux parameter
  //   weibull model:  nvars = 1, ie. shape parameter
  //   gompertz model: nvars = 1, ie. scale parameter
  //   M-spline model: nvars = number of basis terms, ie. spline coefs
  //   B-spline model: nvars = number of basis terms, ie. spline coefs
  vector<lower=0>[Nvars] z_coefs;

  // intercept
  real gamma[N_has_intercept];

  /*
  Not implemented priors
  // parameters for priors
  real<lower=0> global[nhs];
  vector<lower=0>[nhs > 0 ? nK : 0] local[nhs];
  real<lower=0> caux[nhs > 0];
  vector<lower=0>[nK] mix[prior_dist[1] == 5 || prior_dist[1] == 6];
  real<lower=0> ool[prior_dist[1] == 6];

  */
}

transformed parameters {
  // log hazard ratios
  vector[nK] beta;

  // basehaz parameters
  vector[Nvars] coefs;

  // define log hazard ratios
  {
    int pos = 1;
    for (k in 1:nt){
      if(s_K[k] > 0){
        beta[pos:(pos + s_K[k] - 1)] = make_beta(segment(z_beta, pos, s_K[k]), prior_dist[k], segment(prior_mean, pos, s_K[k]), segment(prior_scale, pos, s_K[k]),  segment(prior_df, pos, s_K[k]) ) ;
      }
      pos += s_K[k];
    }
  }


  // define basehaz parameters
  {
    int pos = 1;
    for(k in 1:nt){
      if(s_vars[k] > 0){
        coefs[pos:(pos + s_vars[k] - 1)] = segment(z_coefs, pos, s_vars[k]) .* segment(prior_scale_for_aux, pos, s_vars[k]);
      }
      pos += s_vars[k];
    }
  }
}

model {
  // linear predictor
  vector[Nrcens] eta_rcens;  // time of right censoring
  vector[Nevent] eta_event;  // time of events
 {

  // define iterators
    int pos;
    int pos_e;
    int pos_i_e;
    int pos_rc;
    int pos_i_rc;
    int pos_b;
    int pos_g;
    int pos_coefs;
    int pos_vars;
    int pos_spline_event;
    int pos_spline_rcens;

  // define linear predictor
    pos_e = 1;
    pos_i_e = 1;
    pos_rc = 1;
    pos_i_rc = 1;
    pos_b = 1;
  for(k in 1:nt){
    if(s_K[k] > 0){
      if(s_event[k] > 0){
        matrix[s_event[k], s_K[k]] X_event;
        X_event = to_matrix(
          segment(x_event, pos_e, s_event[k] * s_K[k] ),
          s_event[k], s_K[k] );
          eta_event[pos_i_e:(pos_i_e + s_event[k] - 1) ] = X_event * segment(beta, pos_b, s_K[k]) ;
      }
      if(s_rcens[k] > 0){
        matrix[ s_rcens[k], s_K[k]] X_rcens = to_matrix(
          segment(x_rcens, pos_rc, s_rcens[k] * s_K[k] ),
          s_rcens[k], s_K[k] );
          eta_rcens[pos_i_rc:(pos_i_rc + s_rcens[k] - 1)] = X_rcens * segment(beta, pos_b, s_K[k]);
      }
    } else {
        if (s_event[k] > 0)  eta_event[pos_i_e:(pos_i_e + s_event[k] - 1) ] = rep_vector(0.0, s_event[k]);
        if (s_rcens[k] > 0) eta_rcens[pos_i_rc:(pos_i_rc + s_rcens[k] - 1)] = rep_vector(0.0, s_rcens[k]);
    }
    pos_e += s_event[k] * s_K[k];
    pos_i_e += s_event[k];
    pos_rc += (s_rcens[k] * s_K[k]);
    pos_i_rc += s_rcens[k];
    pos_b += s_K[k];
  }
  // add intercept
  pos_i_e = 1;
  pos_i_rc = 1;
  pos_g = 1;
  for(k in 1:nt){
    if(has_intercept[k] == 1){
      eta_event[pos_i_e : (pos_i_e + s_event[k] - 1)] += gamma[pos_g];
      eta_rcens[pos_i_rc : (pos_i_rc + s_rcens[k] - 1)] += gamma[pos_g];
      pos_g += 1;
    }
    pos_i_e += s_event[k];
    pos_i_rc += s_rcens[k];
  }

  // add on log crude event rate (helps to center intercept)
  pos_i_e = 1;
  pos_i_rc = 1;
  for(k in 1:nt){
     if(s_event[k] > 0) {
       eta_event[pos_i_e:(pos_i_e + s_event[k] - 1)] += log_crude_event_rate[k];
     }
     if(s_rcens[k] > 0) {
       eta_rcens[pos_i_rc:(pos_i_rc + s_rcens[k] - 1)] += log_crude_event_rate[k];
     }
    pos_i_e += s_event[k];
    pos_i_rc += s_rcens[k];
    }
  // evaluate log hazard and log survival
  pos_i_e = 1;
  pos_i_rc = 1;
  pos_coefs = 1;
  pos_spline_event = 1;
  pos_spline_rcens = 1;
  for(k in 1:nt){
    if(type[k] == 1){ // weibull
     real shape = coefs[pos_coefs];
     if(s_event[k] > 0) target += weibull_log_haz( segment(eta_event, pos_i_e, s_event[k]), segment(t_event, pos_i_e, s_event[k]), shape);
     if(s_event[k] > 0) target += weibull_log_surv( segment(eta_event, pos_i_e, s_event[k]), segment(t_event, pos_i_e, s_event[k]), shape);
     if(s_rcens[k] > 0) target += weibull_log_surv( segment(eta_rcens, pos_i_rc, s_rcens[k]), segment(t_rcens, pos_i_rc, s_rcens[k]), shape);

  } else if (type[k] == 4) { // M-splines
      if(s_event[k] > 0){
       matrix[s_event[k], s_vars[k]] iBasis_event = to_matrix(
        segment(ibasis_event, pos_spline_event, s_event[k] * s_vars[k]),
         s_event[k], s_vars[k] );
      matrix[s_event[k], s_vars[k]] Basis_event = to_matrix(
        segment(basis_event, pos_spline_event, s_event[k] * s_vars[k]),
         s_event[k], s_vars[k] );
       target += mspline_log_haz( segment(eta_event, pos_i_e, s_event[k]),
       Basis_event,  segment(coefs, pos_coefs, s_vars[k]) );
       target += mspline_log_surv( segment(eta_event, pos_i_e, s_event[k]), iBasis_event, segment(coefs, pos_coefs, s_vars[k]) ); }
      if(s_rcens[k] > 0){
         matrix[ s_rcens[k], s_vars[k]] iBasis_rcens = to_matrix(
        segment(ibasis_rcens, pos_spline_rcens, s_rcens[k] * s_vars[k]),
         s_rcens[k], s_vars[k] );
         target += mspline_log_surv( segment(eta_rcens, pos_i_rc, s_rcens[k]), iBasis_rcens, segment(coefs, pos_coefs, s_vars[k]) );
        }

      } else if (type[k] == 5) { // exponential model
      if(s_event[k] > 0) target += exponential_log_haz( segment(eta_event, pos_i_e, s_event[k]));
     if(s_event[k] > 0) target += exponential_log_surv( segment(eta_event, pos_i_e, s_event[k]), segment(t_event, pos_i_e, s_event[k]));
     if(s_rcens[k] > 0) target += exponential_log_surv( segment(eta_rcens, pos_i_rc, s_rcens[k]), segment(t_rcens, pos_i_rc, s_rcens[k]));

      } else if (type[k] == 6) { // gompertz model
       real scale = coefs[pos_coefs];
     if(s_event[k] > 0) target += gompertz_log_haz( segment(eta_event, pos_i_e, s_event[k]), segment(t_event, pos_i_e, s_event[k]), scale);
     if(s_event[k] > 0) target += gompertz_log_surv( segment(eta_event, pos_i_e, s_event[k]), segment(t_event, pos_i_e, s_event[k]), scale);
     if(s_rcens[k] > 0) target += gompertz_log_surv( segment(eta_rcens, pos_i_rc, s_rcens[k]), segment(t_rcens, pos_i_rc, s_rcens[k]), scale);
     pos_coefs += 1;

  } else {
    reject("Bug found: invalid baseline hazard (without quadrature).");
  }
  pos_i_e += s_event[k];
  pos_i_rc += s_rcens[k];
  pos_coefs += s_vars[k];
  pos_spline_event += (s_event[k] * s_vars[k]);
  pos_spline_rcens +=  (s_rcens[k] * s_vars[k]);
  }
  //-------- log priors
  // log priors for coefficients
  pos = 1;
  for (k in 1:nt){
    if(s_K[k] > 0){
          real dummy = beta_lp(segment(z_beta, pos, s_K[k]), prior_dist[k], segment(prior_scale, pos, s_K[k]),  segment(prior_df, pos, s_K[k]) ) ;
    }
    pos += s_K[k];
  }
  // log prior for intercept
  pos_g = 1;
  for(k in 1:nt){
    if(has_intercept[k] == 1){
      real dummy = gamma_lp(gamma[pos_g], prior_dist_for_intercept[k],
      prior_mean_for_intercept[k], prior_scale_for_intercept[k], prior_df_for_intercept[k]);
      pos_g += 1;
    }
  }

  // log priors for baseline hazard parameters
  pos = 1;
  for(k in 1:nt){
    if(s_vars[k] > 0){
       real dummy = basehaz_lp(segment(z_coefs, pos, s_vars[k]), prior_dist_for_aux[k], segment(prior_df_for_aux, pos, s_vars[k])); }
  pos += s_vars[k];
  }
 }
}

generated quantities {
  // tranformations to adjust for:
  //   - centering of covariates and
  //   - centering of baseline hazard around the crude event rate
  real alpha[N_has_intercept]; // transformed intercept
  vector[Nvars] aux;              // transformed baseline hazard parameters

  int pos_b = 1;
  int pos_vars = 1;
  int pos_g = 1;
  int pos_int = 1;
  for(k in 1:nt){
    if(type[k] == 4) { // m-splines
      aux[pos_vars:(pos_vars + s_vars[k] - 1)] = segment(coefs, pos_vars, s_vars[k]) * exp(log_crude_event_rate[k] - dot_product(segment(x_bar, pos_b, s_K[k]),
      segment(beta, pos_b, s_K[k])));
    } else { // exp, weibull, gompertz
    aux[pos_vars] = coefs[pos_vars];
    alpha[pos_int] = log_crude_event_rate[k] - dot_product(segment(x_bar, pos_b, s_K[k]), segment(beta, pos_b, s_K[k]) ) + gamma[pos_g];
    pos_g += 1;
    pos_int += 1;
    }
    pos_vars += s_vars[k];
    pos_b += s_K[k];
  }

}

