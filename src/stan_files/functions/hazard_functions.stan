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
