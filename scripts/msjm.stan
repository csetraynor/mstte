functions {
  // Below are code chunks taken from the 'rstanarm' R package,
  // obtained under the terms of the GNU General Public License as
  // published by he Free Software Foundation; either version 3 of
  // the License, or (at your option) any later version.

  /**
   * Create group-specific block-diagonal Cholesky factor, see section 2 of
   * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
   * @param len_theta_L An integer indicating the length of returned vector,
   *   which lme4 denotes as m
   * @param p An integer array with the number variables on the LHS of each |
   * @param dispersion Scalar standard deviation of the errors, calles sigma by lme4
   * @param tau Vector of scale parameters whose squares are proportional to the
   *   traces of the relative covariance matrices of the group-specific terms
   * @param scale Vector of prior scales that are multiplied by elements of tau
   * @param zeta Vector of positive parameters that are normalized into simplexes
   *   and multiplied by the trace of the covariance matrix to produce variances
   * @param rho Vector of radii in the onion method for creating Cholesky factors
   * @param z_T Vector used in the onion method for creating Cholesky factors
   * @return A vector that corresponds to theta in lme4
   */
  vector make_theta_L(int len_theta_L, int[] p, real dispersion,
                      vector tau, vector scale, vector zeta,
                      vector rho, vector z_T) {
    vector[len_theta_L] theta_L;
    int zeta_mark = 1;
    int rho_mark = 1;
    int z_T_mark = 1;
    int theta_L_mark = 1;

    // each of these is a diagonal block of the implicit Cholesky factor
    for (i in 1:size(p)) {
      int nc = p[i];
      if (nc == 1) { // "block" is just a standard deviation
        theta_L[theta_L_mark] = tau[i] * scale[i] * dispersion;
        // unlike lme4, theta[theta_L_mark] includes the dispersion term in it
        theta_L_mark += 1;
      }
      else { // block is lower-triangular
        matrix[nc,nc] T_i;
        real std_dev;
        real T21;
        real trace_T_i = square(tau[i] * scale[i] * dispersion) * nc;
        vector[nc] pi = segment(zeta, zeta_mark, nc); // gamma(zeta | shape, 1)
        pi /= sum(pi);                            // thus dirichlet(pi | shape)

        // unlike lme4, T_i includes the dispersion term in it
        zeta_mark += nc;
        std_dev = sqrt(pi[1] * trace_T_i);
        T_i[1,1] = std_dev;

        // Put a correlation into T_i[2,1] and scale by std_dev
        std_dev = sqrt(pi[2] * trace_T_i);
        T21 = 2.0 * rho[rho_mark] - 1.0;
        rho_mark += 1;
        T_i[2,2] = std_dev * sqrt(1.0 - square(T21));
        T_i[2,1] = std_dev * T21;

        for (r in 2:(nc - 1)) { // scaled onion method to fill T_i
          int rp1 = r + 1;
          vector[r] T_row = segment(z_T, z_T_mark, r);
          real scale_factor = sqrt(rho[rho_mark] / dot_self(T_row)) * std_dev;
          z_T_mark += r;
          std_dev = sqrt(pi[rp1] * trace_T_i);
          for(c in 1:r) T_i[rp1,c] = T_row[c] * scale_factor;
          T_i[rp1,rp1] = sqrt(1.0 - rho[rho_mark]) * std_dev;
          rho_mark += 1;
        }

        // now vech T_i
        for (c in 1:nc) for (r in c:nc) {
          theta_L[theta_L_mark] = T_i[r,c];
          theta_L_mark += 1;
        }
      }
    }
    return theta_L;
  }

  /**
  * Create group-specific coefficients, see section 2 of
  * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  *
  * @param z_b Vector whose elements are iid normal(0,sigma) a priori
  * @param theta Vector with covariance parameters as defined in lme4
  * @param p An integer array with the number variables on the LHS of each |
  * @param l An integer array with the number of levels for the factor(s) on
  *   the RHS of each |
  * @return A vector of group-specific coefficients
  */
  vector make_b(vector z_b, vector theta_L, int[] p, int[] l) {
    vector[rows(z_b)] b;
    int b_mark = 1;
    int theta_L_mark = 1;
    for (i in 1:size(p)) {
      int nc = p[i];
      if (nc == 1) {
        real theta_L_start = theta_L[theta_L_mark];
        for (s in b_mark:(b_mark + l[i] - 1))
          b[s] = theta_L_start * z_b[s];
        b_mark += l[i];
        theta_L_mark += 1;
      }
      else {
        matrix[nc,nc] T_i = rep_matrix(0, nc, nc);
        for (c in 1:nc) {
          T_i[c,c] = theta_L[theta_L_mark];
          theta_L_mark += 1;
          for(r in (c+1):nc) {
            T_i[r,c] = theta_L[theta_L_mark];
            theta_L_mark += 1;
          }
        }
        for (j in 1:l[i]) {
          vector[nc] temp = T_i * segment(z_b, b_mark, nc);
          b_mark -= 1;
          for (s in 1:nc) b[b_mark + s] = temp[s];
          b_mark += nc + 1;
        }
      }
    }
    return b;
  }

  /**
   * Prior on group-specific parameters
   *
   * @param z_b A vector of primitive coefficients
   * @param z_T A vector of primitives for the unit vectors in the onion method
   * @param rho A vector radii for the onion method
   * @param zeta A vector of primitives for the simplexes
   * @param tau A vector of scale parameters
   * @param regularization A real array of LKJ hyperparameters
   * @param delta A real array of concentration paramters
   * @param shape A vector of shape parameters
   * @param t An integer indicating the number of group-specific terms
   * @param p An integer array with the number variables on the LHS of each |
   * @return target()
   */
  real decov_lp(vector z_b, vector z_T, vector rho, vector zeta, vector tau,
                real[] regularization, real[] delta, vector shape,
                int t, int[] p) {
    int pos_reg = 1;
    int pos_rho = 1;
    target += normal_lpdf(z_b | 0, 1);
    target += normal_lpdf(z_T | 0, 1);
    for (i in 1:t) if (p[i] > 1) {
      vector[p[i] - 1] shape1;
      vector[p[i] - 1] shape2;
      real nu = regularization[pos_reg] + 0.5 * (p[i] - 2);
      pos_reg += 1;
      shape1[1] = nu;
      shape2[1] = nu;
      for (j in 2:(p[i]-1)) {
        nu -= 0.5;
        shape1[j] = 0.5 * j;
        shape2[j] = nu;
      }
      target += beta_lpdf(rho[pos_rho:(pos_rho + p[i] - 2)] | shape1, shape2);
      pos_rho += p[i] - 1;
    }
    target += gamma_lpdf(zeta | delta, 1);
    target += gamma_lpdf(tau  | shape, 1);
    return target();
  }

  /**
   * Hierarchical shrinkage parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @param global_prior_scale A positive real number
   * @param error_scale 1 or sigma in the Gaussian case
   * @param c2 A positive real number
   * @return A vector of coefficientes
   */
  vector hs_prior(vector z_beta, real[] global, vector[] local,
                  real global_prior_scale, real error_scale, real c2) {
    int K = rows(z_beta);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt( c2 * lambda2 ./ (c2 + square(tau) * lambda2) );
    return z_beta .* lambda_tilde * tau;
  }

  /**
   * Hierarchical shrinkage plus parameterization
   *
   * @param z_beta A vector of primitive coefficients
   * @param global A real array of positive numbers
   * @param local A vector array of positive numbers
   * @param global_prior_scale A positive real number
   * @param error_scale 1 or sigma in the Gaussian case
   * @param c2 A positive real number
   * @return A vector of coefficientes
   */
  vector hsplus_prior(vector z_beta, real[] global, vector[] local,
                      real global_prior_scale, real error_scale, real c2) {
    int K = rows(z_beta);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    vector[K] eta = local[3] .* sqrt(local[4]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda_eta2 = square(lambda .* eta);
    vector[K] lambda_tilde = sqrt( c2 * lambda_eta2 ./
                                 ( c2 + square(tau) * lambda_eta2) );
    return z_beta .* lambda_tilde * tau;
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
   * Return two-dimensional array of group membership
   *
   * @param N An integer indicating the number of observations
   * @param t An integer indicating the number of grouping variables
   * @param v An integer array with the indices of group membership
   * @return An two-dimensional integer array of group membership
   */
  int[,] make_V(int N, int t, int[] v) {
    int V[t,N];
    int pos = 1;
    if (t > 0) for (j in 1:N) for (i in 1:t) {
      V[i,j] = v[pos] + 1;
      pos += 1;
    }
    return V;
  }

  /**
  * faster version of csr_matrix_times_vector
  * declared here and defined in C++
  *
  * @param m Integer number of rows
  * @param n Integer number of columns
  * @param w Vector (see reference manual)
  * @param v Integer array (see reference manual)
  * @param u Integer array (see reference manual)
  * @param b Vector that is multiplied from the left by the CSR matrix
  * @return A vector that is the product of the CSR matrix and b
  */
  vector csr_matrix_times_vector2(int m, int n, vector w,
                                  int[] v, int[] u, vector b);

  /**
   * Calculate lower bound on intercept
   *
   * @param family Integer family code
   *   1 = gaussian
   *   2 = gamma
   *   3 = inv-gaussian
   *   4 = beta
   *   5 = binomial
   *   6 = poisson
   *   7 = neg-binom
   *   8 = poisson w/ gamma noise (not currently used but in count.stan)
   * @param link Integer link code
   * @return real lower bound
   */
  real make_lower(int family, int link) {
    if (family == 1) return negative_infinity(); // Gaussian
    if (family <= 3) { // Gamma or inverse Gaussian
      if (link == 2) return negative_infinity(); // log
      return 0;
    }
    return negative_infinity();
  }

  /**
   * Calculate upper bound on intercept
   *
   * @param family Integer family code (see make_lower above for codes)
   * @param link Integer link code
   * @return real upper bound
   */
  real make_upper(int family, int link) {
    if (family == 4 && link == 5) return 0;
    return positive_infinity();
  }

    /**
   * Apply inverse link function to linear predictor
   * see help(binom) in R
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_bern(vector eta, int link) {
    if (link == 1)      return(inv_logit(eta)); // logit
    else if (link == 2) return(Phi(eta)); // probit
    else if (link == 3) return(atan(eta) / pi() + 0.5);  // cauchit
    else if (link == 4) return(exp(eta)); // log
    else if (link == 5) return(inv_cloglog(eta)); // cloglog
    else reject("Invalid link");
    return eta; // never reached
  }

  /**
   * Increment with the unweighted log-likelihood
   * @param link An integer indicating the link function
   * @param eta0 A vector of linear predictors | y = 0
   * @param eta1 A vector of linear predictors | y = 1
   * @param N An integer array of length 2 giving the number of
   *   observations where y = 0 and y = 1 respectively
   * @return lp__
   */
  real ll_bern_lp(vector eta0, vector eta1, int link, int[] N) {
    if (link == 1) { // logit
      target += logistic_lccdf(eta0 | 0, 1);
      target += logistic_lcdf( eta1 | 0, 1);
    }
    else if (link == 2) {  // probit
      target += normal_lccdf(eta0 | 0, 1);
      target += normal_lcdf( eta1 | 0, 1);
    }
    else if (link == 3) {  // cauchit
      target += cauchy_lccdf(eta0 | 0, 1);
      target += cauchy_lcdf( eta1 | 0, 1);
    }
    else if(link == 4) {  // log
      target += log1m_exp(eta0);
      target += eta1;  // already in log form
    }
    else if(link == 5) {  // cloglog
      target += log1m_exp(-exp(eta1));
      target += -exp(eta0);
    }
    else reject("Invalid link");
    return target();
  }

  /**
   * Pointwise (pw) log-likelihood vector
   *
   * @param y The integer outcome variable. Note that function is
   *  called separately with y = 0 and y = 1
   * @param eta Vector of linear predictions
   * @param link An integer indicating the link function
   * @return A vector
   */
  vector pw_bern(int y, vector eta, int link) {
    int N = rows(eta);
    vector[N] ll;
    if (link == 1) {  // logit
      for (n in 1:N) ll[n] = bernoulli_logit_lpmf(y | eta[n]);
    }
    else if (link <= 5) {  // link = probit, cauchit, log, or cloglog
      vector[N] pi = linkinv_bern(eta, link); // may not be stable
      for (n in 1:N) ll[n] = bernoulli_lpmf(y | pi[n]);
    }
    else reject("Invalid link");
    return ll;
  }

  /**
   * Log-normalizing constant in the clogit case
   *
   * @param N_j Integer number of observations in the j-th group
   * @param D_j Integer number of successes in the j-th group
   * @param eta_j Vector of linear predictions in the j-th group
   * @return A scalar that normalizes the probabilities on the log-scale
   */
  real log_clogit_denom(int N_j, int D_j, vector eta_j);
  real log_clogit_denom(int N_j, int D_j, vector eta_j) {
    if (D_j == 1 && N_j == rows(eta_j)) return log_sum_exp(eta_j);
    if (D_j == 0) return 0;
    if (N_j == D_j) {
      if (D_j == 1) return eta_j[N_j];
      return sum(segment(eta_j, N_j - 1, 2));
    }
    else {
      int N_jm1 = N_j - 1;
      return log_sum_exp(log_clogit_denom(N_jm1, D_j, eta_j),
                         log_clogit_denom(N_jm1, D_j - 1, eta_j) + eta_j[N_j]);
    }
    return not_a_number();  // never reaches
  }

  /**
   * Log-likelihood for a clogit model
   * @param eta0 Linear predictors when y == 0
   * @param eta1 Linear predictors when y == 1
   * @param successes Integer array with the number of successes in group j
   * @param failures Integer array with the number of failures in group j
   * @param observations Integer array with the number of observations in group j
   * @return lp__
   */
  real ll_clogit_lp(vector eta0, vector eta1,
                    int[] successes, int[] failures, int[] observations) {
    int J = num_elements(observations);
    int pos0 = 1;
    int pos1 = 1;
    vector[J] summands;
    for (j in 1:J) {
      int D_g = successes[j];
      int N_g = observations[j];
      int F_g = failures[j];
      vector[N_g] eta_g = append_row(segment(eta1, pos1, D_g),
                                     segment(eta0, pos0, F_g));
      summands[j] = log_clogit_denom(N_g, D_g, eta_g);
      pos0 += F_g;
      pos1 += D_g;
    }
    target += sum(eta1) - sum(summands);
    return target();
  }
  /**
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_binom(vector eta, int link) {
    if (link == 1)      return(inv_logit(eta)); // logit
    else if (link == 2) return(Phi(eta)); // probit
    else if (link == 3) return(atan(eta) / pi() + 0.5);  // cauchit
    else if (link == 4) return(exp(eta)); // log
    else if (link == 5) return(inv_cloglog(eta)); // cloglog
    else reject("Invalid link");
    return eta; // never reached
  }

  /**
  * Increment with the unweighted log-likelihood
  * @param y An integer array indicating the number of successes
  * @param trials An integer array indicating the number of trials
  * @param eta A vector of linear predictors
  * @param link An integer indicating the link function
  * @return lp__
  */
  real ll_binom_lp(int[] y, int[] trials, vector eta, int link) {
    if (link == 1) target += binomial_logit_lpmf(y | trials, eta);
    else if (link <  4) target += binomial_lpmf( y | trials, linkinv_binom(eta, link));
    else if (link == 4) {  // log
      for (n in 1:num_elements(y)) {
        target += y[n] * eta[n];
        target += (trials[n] - y[n]) * log1m_exp(eta[n]);
        target += lchoose(trials[n], y[n]);
      }
    }
    else if (link == 5) {  // cloglog
      for (n in 1:num_elements(y)) {
        real neg_exp_eta = -exp(eta[n]);
        target += y[n] * log1m_exp(neg_exp_eta);
        target += (trials[n] - y[n]) * neg_exp_eta;
        target += lchoose(trials[n], y[n]);
      }
    }
    else reject("Invalid link");
    return target();
  }

  /**
  * Pointwise (pw) log-likelihood vector
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_binom(int[] y, int[] trials, vector eta, int link) {
    int N = rows(eta);
    vector[N] ll;
    if (link == 1) {  // logit
      for (n in 1:N)
        ll[n] = binomial_logit_lpmf(y[n] | trials[n], eta[n]);
    }
    else if (link <= 5) {  // link = probit, cauchit, log, or cloglog
      vector[N] pi = linkinv_binom(eta, link); // may be unstable
      for (n in 1:N) ll[n] = binomial_lpmf(y[n] | trials[n], pi[n]) ;
    }
    else reject("Invalid link");
    return ll;
  }
  /**
   * Apply inverse link function to linear predictor
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_gauss(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta);
    else if (link == 3) return inv(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /**
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gauss(vector y, vector eta, real sigma, int link) {
    return -0.5 * log(6.283185307179586232 * sigma) -
            0.5 * square((y - linkinv_gauss(eta, link)) / sigma);
  }

  /**
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_gamma(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta);
    else if (link == 3) return inv(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /**
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param eta A vector of linear predictors
  * @param shape A real number for the shape parameter
  * @param link An integer indicating the link function
  * @param sum_log_y A scalar equal to the sum of log(y)
  * @return A scalar log-likelihood
  */
  real GammaReg(vector y, vector eta, real shape,
                int link, real sum_log_y) {
    real ret = rows(y) * (shape * log(shape) - lgamma(shape)) +
               (shape - 1) * sum_log_y;
    if (link == 2)      // link is log
      ret -= shape * sum(eta) + shape * sum(y ./ exp(eta));
    else if (link == 1) // link is identity
      ret -= shape * sum(log(eta)) + shape * sum(y ./ eta);
    else if (link == 3) // link is inverse
      ret += shape * sum(log(eta)) - shape * dot_product(eta, y);
    else reject("Invalid link");
    return ret;
  }

  /**
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param shape A real number for the shape parameter
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_gamma(vector y, vector eta, real shape, int link) {
    int N = rows(eta);
    vector[N] ll;
    if (link == 3) { // link = inverse
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape * eta[n]);
      }
    }
    else if (link == 2) { // link = log
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape / exp(eta[n]));
      }
    }
    else if (link == 1) { // link = identity
      for (n in 1:N) {
        ll[n] = gamma_lpdf(y[n] | shape, shape / eta[n]);
      }
    }
    else reject("Invalid link");
    return ll;
  }

  /**
  * Apply inverse link function to linear predictor
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_inv_gaussian(vector eta, int link) {
    if (link == 1)      return eta;
    else if (link == 2) return exp(eta);
    else if (link == 3) return inv(eta);
    else if (link == 4) return inv_sqrt(eta);
    else reject("Invalid link");
    return eta; // never reached
  }

  /**
  * inverse Gaussian log-PDF
  *
  * @param y The vector of outcomes
  * @param mu The vector of conditional means
  * @param lambda A positive scalar dispersion parameter
  * @param sum_log_y A scalar equal to the sum of log(y)
  * @param sqrt_y A vector equal to sqrt(y)
  * @return A scalar
  */
  real inv_gaussian(vector y, vector mu, real lambda,
                    real sum_log_y, vector sqrt_y) {
    return 0.5 * rows(y) * log(lambda / 6.283185307179586232) -
      1.5 * sum_log_y -
      0.5 * lambda * dot_self( (y - mu) ./ (mu .* sqrt_y) );
  }

  /**
  * Pointwise (pw) log-likelihood vector
  *
  * @param y A vector corresponding to the outcome variable.
  * @param eta The linear predictors
  * @param lamba A positive scalar dispersion parameter
  * @param link An integer indicating the link function
  * @param log_y A precalculated vector of the log of y
  * @param sqrt_y A precalculated vector of the square root of y
  * @return A vector of log-likelihoods
  */
  vector pw_inv_gaussian(vector y, vector eta, real lambda,
                         int link, vector log_y, vector sqrt_y) {
    vector[rows(y)] mu = linkinv_inv_gaussian(eta, link); // link checked
    return -0.5 * lambda * square( (y - mu) ./ (mu .* sqrt_y) ) +
            0.5 * log(lambda / 6.283185307179586232) - 1.5 * log_y;
  }

  /**
  * PRNG for the inverse Gaussian distribution
  *
  * Algorithm from wikipedia
  *
  * @param mu The expectation
  * @param lambda The dispersion
  * @return A draw from the inverse Gaussian distribution
  */
  real inv_gaussian_rng(real mu, real lambda) {
    real mu2 = square(mu);
    real z = uniform_rng(0,1);
    real y = square(normal_rng(0,1));
    real x = mu + ( mu2 * y - mu * sqrt(4 * mu * lambda * y + mu2 * square(y)) )
           / (2 * lambda);
    if (z <= (mu / (mu + x))) return x;
    else return mu2 / x;
  }

  /**
  * Apply inverse link function to linear predictor for beta models
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_beta(vector eta, int link) {
    if (link == 1) return inv_logit(eta);  // logit
    else if (link == 2) return Phi(eta);   // probit
    else if (link == 3) return inv_cloglog(eta);  // cloglog
    else if (link == 4) return 0.5 + atan(eta) / pi(); // cauchy
    else if (link == 5) return exp(eta); // log
    else if (link == 6) return 1 - inv_cloglog(-eta); // loglog
    else reject("invalid link");
    return eta; // never reached
  }

  /**
  * Apply inverse link function to linear predictor for dispersion for beta models
  *
  * @param eta Linear predictor vector
  * @param link An integer indicating the link function
  * @return A vector, i.e. inverse-link(eta)
  */
  vector linkinv_beta_z(vector eta, int link) {
    if (link == 1) return exp(eta);         // log
    else if (link == 2) return eta;         // identity
    else if (link == 3) return square(eta); // sqrt
    else reject("Invalid link")
    return eta; // never reached
  }

  /**
  * Pointwise (pw) log-likelihood vector for beta models
  *
  * @param y The vector of outcomes
  * @param eta The linear predictors
  * @param dispersion Positive dispersion parameter
  * @param link An integer indicating the link function
  * @return A vector of log-likelihoods
  */
  vector pw_beta(vector y, vector eta, real dispersion, int link) {
    vector[rows(y)] ll;
    vector[rows(y)] mu = linkinv_beta(eta, link); // link checked
    for (n in 1:rows(y)) {
      ll[n] = beta_lpdf(y[n] | mu[n] * dispersion, (1 - mu[n]) * dispersion);
    }
    return ll;
  }

  /**
  * Pointwise (pw) log-likelihood vector for beta models with z variables
  *
  * @param y The vector of outcomes
  * @param eta The linear predictors (for y)
  * @param eta_z The linear predictors (for dispersion)
  * @param link An integer indicating the link function passed to linkinv_beta
  * @param link_phi An integer indicating the link function passed to linkinv_beta_z
  * @return A vector of log-likelihoods
  */
  vector pw_beta_z(vector y, vector eta, vector eta_z, int link, int link_phi) {
    vector[rows(y)] ll;
    vector[rows(y)] mu = linkinv_beta(eta, link); // link checked
    vector[rows(y)] mu_z = linkinv_beta_z(eta_z, link_phi); // link checked
    for (n in 1:rows(y)) {
      ll[n] = beta_lpdf(y[n] | mu[n] * mu_z[n], (1-mu[n]) * mu_z[n]);
    }
    return ll;
  }
  /**
   * Apply inverse link function to linear predictor
   * see help(poisson) in R
   *
   * @param eta Linear predictor vector
   * @param link An integer indicating the link function
   * @return A vector, i.e. inverse-link(eta)
   */
  vector linkinv_count(vector eta, int link) {
    if (link == 1)      return exp(eta);     // log
    else if (link == 2) return eta;          // identity
    else if (link == 3) return(square(eta)); // sqrt
    else reject("Invalid link");
    return eta; // never reached
  }

  /**
  * Pointwise (pw) log-likelihood vector for the Poisson distribution
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param eta The vector of linear predictors
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_pois(int[] y, vector eta, int link) {
    int N = rows(eta);
    vector[N] ll;
    if (link == 1)  // log
      for (n in 1:N) ll[n] = poisson_log_lpmf(y[n] | eta[n]);
    else if (link <= 3) {  // link = identity or sqrt
      vector[N] phi = linkinv_count(eta, link);
      for (n in 1:N) ll[n] = poisson_lpmf(y[n] | phi[n]) ;
    }
    else reject("Invalid link");
    return ll;
  }

  /**
  * Pointwise (pw) log-likelihood vector for the negative binomial distribution
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param eta The vector of linear predictors
  * @param theta The reciprocal_dispersion parameter
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_nb(int[] y, vector eta, real theta, int link) {
    int N = rows(eta);
    vector[N] rho = linkinv_count(eta, link); // link checked
    vector[N] ll;
    for (n in 1:N) ll[n] = neg_binomial_2_lpmf(y[n] | rho[n], theta);
    return ll;
  }
  /**
  * Return the required number of local hs parameters
  *
  * @param prior_dist An integer indicating the prior distribution
  * @return An integer
  */
  int get_nvars_for_hs(int prior_dist) {
    int hs = 0;
    if (prior_dist == 3) hs = 2;
    else if (prior_dist == 4) hs = 4;
    return hs;
  }

  /**
  * Return the lower/upper bound for the specified intercept type
  *
  * @param intercept_type An integer specifying the type of intercept;
  *   0=no intercept, 1=unbounded, 2=lower bounded, 3=upper bounded
  * @return A real, corresponding to the lower bound
  */
  real lb(int intercept_type) {
    real lb;
    if (intercept_type == 2) lb = 0;
    else lb = negative_infinity();
    return lb;
  }
  real ub(int intercept_type) {
    real ub;
    if (intercept_type == 3) ub = 0;
    else ub = positive_infinity();
    return ub;
  }

  /**
  * Get the indices corresponding to the lower tri of a square matrix
  *
  * @param dim The number of rows in the square matrix
  * @return A vector of indices
  */
  int[] lower_tri_indices(int dim) {
    int indices[dim + choose(dim, 2)];
    int mark = 1;
    for (r in 1:dim) {
      for (c in r:dim) {
        indices[mark] = (r - 1) * dim + c;
        mark += 1;
      }
    }
    return indices;
  }

  /**
  * Scale the auxiliary parameter based on prior information
  *
  * @param aux_unscaled A real, the unscaled auxiliary parameter
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Real scalars, the mean and scale
  *   of the prior distribution
  * @return A real, corresponding to the scaled auxiliary parameter
  */
  real make_aux(real aux_unscaled, int prior_dist,
                real prior_mean, real prior_scale) {
    real aux;
    if (prior_dist == 0) // none
      aux = aux_unscaled;
    else {
      aux = prior_scale * aux_unscaled;
      if (prior_dist <= 2) // normal or student_t
        aux += prior_mean;
    }
    return aux;
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
                   vector prior_scale, vector prior_df, real global_prior_scale,
                   real[] global, vector[] local, real[] ool, vector[] mix,
                   real[] aux, int family, real slab_scale, real[] caux) {
    vector[rows(z_beta)] beta;
    if (prior_dist == 0) beta = z_beta;
    else if (prior_dist == 1) beta = z_beta .* prior_scale + prior_mean;
    else if (prior_dist == 2) for (k in 1:rows(prior_mean)) {
      beta[k] = CFt(z_beta[k], prior_df[k]) * prior_scale[k] + prior_mean[k];
    }
    else if (prior_dist == 3) {
      real c2 = square(slab_scale) * caux[1];
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hs_prior(z_beta, global, local, global_prior_scale, aux[1], c2);
      else
        beta = hs_prior(z_beta, global, local, global_prior_scale, 1, c2);
    }
    else if (prior_dist == 4) {
      real c2 = square(slab_scale) * caux[1];
      if (family == 1) // don't need is_continuous since family == 1 is gaussian in mvmer
        beta = hsplus_prior(z_beta, global, local, global_prior_scale, aux[1], c2);
      else
        beta = hsplus_prior(z_beta, global, local, global_prior_scale, 1, c2);
    }
    else if (prior_dist == 5) // laplace
      beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    else if (prior_dist == 6) // lasso
      beta = prior_mean + ool[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    return beta;
  }

  /**
  * Create group-specific coefficients, see section 2 of
  * https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  *
  * @param z_b Vector whose elements are iid normal(0,sigma) a priori
  * @param theta Vector with covariance parameters as defined in lme4
  * @param p An integer array with the number variables on the LHS of each |
  * @param l An integer array with the number of levels for the factor(s) on
  *   the RHS of each |
  * @param i The index of the grouping factor for which you want to return
  *   the group-specific coefficients for
  * @return An array of group-specific coefficients for grouping factor i
  */
  matrix make_b_matrix(vector z_b, vector theta_L, int[] p, int[] l, int i) {
    matrix[p[i],l[i]] b_matrix;
    int nc = p[i];
    int b_mark = 1;
    int theta_L_mark = 1;
    if (i > 1) {
      for (j in 1:(i-1)) {
        theta_L_mark += p[j] + choose(p[j], 2);
        b_mark += p[j] * l[j];
      }
    }
    if (nc == 1) {
      real theta_L_start = theta_L[theta_L_mark];
      for (s in b_mark:(b_mark + l[i] - 1))
        b_matrix[nc,s] = theta_L_start * z_b[s];
    }
    else {
      matrix[nc,nc] T_i = rep_matrix(0, nc, nc);
      for (c in 1:nc) {
        T_i[c,c] = theta_L[theta_L_mark];
        theta_L_mark += 1;
        for(r in (c+1):nc) {
          T_i[r,c] = theta_L[theta_L_mark];
          theta_L_mark += 1;
        }
      }
      for (j in 1:l[i]) {
        vector[nc] temp = T_i * segment(z_b, b_mark, nc);
        b_matrix[,j] = temp;
        b_mark += nc;
      }
    }
    return b_matrix';
  }

  /**
  * Evaluate the linear predictor for the glmer submodel
  *
  * @param X Design matrix for fe
  * @param Z1 Design matrix for re, for first grouping factor
  * @param Z2 Design matrix for re, for second grouping factor
  * @param Z1_id Group indexing for Z1
  * @param Z2_id Group indexing for Z2
  * @param gamma The intercept parameter
  * @param beta Vector of population level parameters
  * @param b1Mat Matrix of group level params for first grouping factor
  * @param b2Mat Matrix of group level params for second grouping factor
  * @param b1Mat_colshift,b2Mat_colshift Number of columns in b1Mat/b2Mat
  *   that correpond to group level params from prior glmer submodels
  * @param intercept_type The type of intercept parameter (0 = none,
  *   1 = unbounded, 2 = lower bound, 3 = upper bound)
  * @return A vector containing the linear predictor for the glmer submodel
  */
  vector evaluate_eta(matrix X, vector[] Z1, vector[] Z2, int[] Z1_id, int[] Z2_id,
                      real[] gamma, vector beta, matrix b1Mat, matrix b2Mat,
                      int b1Mat_colshift, int b2Mat_colshift,
                      int intercept_type) {
    int N = rows(X);    // num rows in design matrix
    int K = rows(beta); // num predictors
    int p1 = size(Z1);  // num group level params for group factor 1
    int p2 = size(Z2);  // num group level params for group factor 2
    vector[N] eta;

    if (K > 0) eta = X * beta;
    else eta = rep_vector(0.0, N);

    if (intercept_type > 0) { // submodel has an intercept
      if (intercept_type == 1) eta += gamma[1];
      else if (intercept_type == 2) eta += gamma[1] - max(eta);
      else if (intercept_type == 3) eta += gamma[1] - min(eta);
    }

    if (p1 > 0) { // submodel includes group factor 1
      for (k in 1:p1)
        for (n in 1:N)
          eta[n] += (b1Mat[Z1_id[n], k+b1Mat_colshift]) * Z1[k,n];
    }
    if (p2 > 0) { // submodel includes group factor 2
      for (k in 1:p2)
        for (n in 1:N)
          eta[n] += (b2Mat[Z2_id[n], k+b2Mat_colshift]) * Z2[k,n];
    }

    return eta;
  }

  /**
  * Evaluate mu based on eta, family and link
  *
  * @param eta Vector of linear predictors
  * @param family An integer indicating the family
  * @param link An integer indicating the link function (differs by family)
  * @return A vector
  */
  vector evaluate_mu(vector eta, int family, int link) {
    vector[rows(eta)] mu;
    if (family == 1)
      mu = linkinv_gauss(eta, link);
    else if (family == 2)
      mu = linkinv_gamma(eta, link);
    else if (family == 3)
      mu = linkinv_inv_gaussian(eta, link);
    else if (family == 4)
      mu = linkinv_bern(eta, link);
    else if (family == 5)
      mu = linkinv_binom(eta, link);
    else if (family == 6 || family == 7 || family == 8)
      mu = linkinv_count(eta, link);
    return mu;
  }

  /**
  * Increment the target with the log-likelihood for the glmer submodel
  *
  * @param z_beta A vector of primitive parameters
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Vectors of mean and scale parameters
  *   for the prior distributions
  * @return A vector containing the population level parameters (coefficients)
  */
  void glm_lp(vector y_real, int[] y_integer, vector eta, real[] aux,
              int family, int link, real sum_log_y, vector sqrt_y, vector log_y) {
    if (family == 1) {  // gaussian
      if (link == 1) target += normal_lpdf(y_real | eta, aux[1]);
      else if (link == 2) target += lognormal_lpdf(y_real | eta, aux[1]);
      else target += normal_lpdf(y_real | inv(eta), aux[1]);
    }
    else if (family == 2) {  // gamma
      target += GammaReg(y_real, eta, aux[1], link, sum_log_y);
    }
    else if (family == 3) {  // inverse gaussian
      target += inv_gaussian(y_real, linkinv_inv_gaussian(eta, link),
                             aux[1], sum_log_y, sqrt_y);
    }
    else if (family == 4) {  // bernoulli
      if (link == 1) target += bernoulli_logit_lpmf(y_integer | eta);
      else target += bernoulli_lpmf(y_integer | linkinv_bern(eta, link));
    }
    else if (family == 5) {  // binomial
      reject("Binomial with >1 trials not allowed.");
    }
    else if (family == 6 || family == 8) {  // poisson or poisson-gamma
      if (link == 1) target += poisson_log_lpmf(y_integer | eta);
      else target += poisson_lpmf(y_integer | linkinv_count(eta, link));
    }
    else if (family == 7) {  // negative binomial
        if (link == 1) target += neg_binomial_2_log_lpmf(y_integer | eta, aux[1]);
      else target += neg_binomial_2_lpmf(y_integer | linkinv_count(eta, link), aux[1]);
    }
    else reject("Invalid family.");
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
  * @return nothing
  */
  void beta_lp(vector z_beta, int prior_dist, vector prior_scale,
               vector prior_df, real global_prior_df, vector[] local,
               real[] global, vector[] mix, real[] one_over_lambda,
               real slab_df, real[] caux) {
    if      (prior_dist == 1) target += normal_lpdf(z_beta | 0, 1);
    else if (prior_dist == 2) target += normal_lpdf(z_beta | 0, 1); // Student t
    else if (prior_dist == 3) { // hs
      target += normal_lpdf(z_beta | 0, 1);
      target += normal_lpdf(local[1] | 0, 1);
      target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      target += normal_lpdf(global[1] | 0, 1);
      target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
      target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
    }
    else if (prior_dist == 4) { // hs+
      target += normal_lpdf(z_beta | 0, 1);
      target += normal_lpdf(local[1] | 0, 1);
      target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
      target += normal_lpdf(local[3] | 0, 1);
      // unorthodox useage of prior_scale as another df hyperparameter
      target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);
      target += normal_lpdf(global[1] | 0, 1);
      target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
      target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);
    }
    else if (prior_dist == 5) { // laplace
      target += normal_lpdf(z_beta | 0, 1);
      target += exponential_lpdf(mix[1] | 1);
    }
    else if (prior_dist == 6) { // lasso
      target += normal_lpdf(z_beta | 0, 1);
      target += exponential_lpdf(mix[1] | 1);
      target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);
    }
    else if (prior_dist == 7) { // product_normal
      target += normal_lpdf(z_beta | 0, 1);
    }
    /* else prior_dist is 0 and nothing is added */
  }

  /**
  * Log-prior for intercept parameters
  *
  * @param gamma Real, the intercept parameter
  * @param dist Integer, the type of prior distribution
  * @param mean Real, mean of prior distribution
  * @param scale Real, scale for the prior distribution
  * @param df Real, df for the prior distribution
  * @return nothing
  */
  void gamma_lp(real gamma, int dist, real mean, real scale, real df) {
    if (dist == 1)  // normal
      target += normal_lpdf(gamma | mean, scale);
    else if (dist == 2)  // student_t
      target += student_t_lpdf(gamma | df, mean, scale);
    /* else dist is 0 and nothing is added */
  }

  /**
  * Log-prior for auxiliary parameters
  *
  * @param aux_unscaled Vector (potentially of length 1) of unscaled
  *   auxiliary parameter(s)
  * @param dist Integer specifying the type of prior distribution
  * @param scale Real specifying the scale for the prior distribution
  * @param df Real specifying the df for the prior distribution
  * @return nothing
  */
  void aux_lp(real aux_unscaled, int dist, real scale, real df) {
    if (dist > 0 && scale > 0) {
      if (dist == 1)
        target += normal_lpdf(aux_unscaled | 0, 1);
      else if (dist == 2)
        target += student_t_lpdf(aux_unscaled | df, 0, 1);
      else
        target += exponential_lpdf(aux_unscaled | 1);
    }
  }

  /**
  * Evaluate the mean of the posterior predictive distribution
  *
  * @param mu Vector containing the mean of the posterior predictive
  *   distribution for each observation (ie. the linear predictor after
  *   applying the inverse link function).
  * @param real The auxiliary parameter for the glmer submodel. This will be
  *   an empty array if the submodel does not have an auxiliary parameter
  * @param family An integer specifying the family
  * @return A real, the mean of the posterior predictive distribution
  */
  real mean_PPD_rng(vector mu, real[] aux, int family) {
    int N = rows(mu);
    real mean_PPD = 0;
    if (family == 1) { // gaussian
      for (n in 1:N)
        mean_PPD += normal_rng(mu[n], aux[1]);
    }
    else if (family == 2) {  // gamma
      for (n in 1:N)
        mean_PPD += gamma_rng(aux[1], aux[1] / mu[n]);
    }
    else if (family == 3) {  // inverse gaussian
      for (n in 1:N)
        mean_PPD += inv_gaussian_rng(mu[n], aux[1]);
    }
    else if (family == 4) {  // bernoulli
      for (n in 1:N)
        mean_PPD += bernoulli_rng(mu[n]);
    }
    else if (family == 5) {  // binomial
      reject("Binomial with >1 trials not allowed.");
    }
    else if (family == 6 || family == 8) {
      real poisson_max = pow(2.0, 30.0);
      for (n in 1:N) {  // poisson or poisson-gamma
        if (mu[n] < poisson_max)
          mean_PPD += poisson_rng(mu[n]);
        else
          mean_PPD += normal_rng(mu[n], sqrt(mu[n]));
      }
    }
    else if (family == 7) {
      real poisson_max = pow(2.0, 30.0);
      for (n in 1:N) {  // negative binomial
        real gamma_temp;
        if (is_inf(aux[1]))
          gamma_temp = mu[n];
        else
          gamma_temp = gamma_rng(aux[1], aux[1] / mu[n]);
        if (gamma_temp < poisson_max)
          mean_PPD += poisson_rng(gamma_temp);
        else
          mean_PPD += normal_rng(gamma_temp, sqrt(gamma_temp));
      }
    }
    return mean_PPD / N;
  }

    /**
  * Return the lower bound for the baseline hazard parameters
  *
  * @param type An integer indicating the type of baseline hazard
  * @return A real
  */
  real coefs_lb(int type) {
    real lbound;
    if (type == 2) // B-splines, on log haz scale
      lbound = negative_infinity();
    else if (type == 3) // piecewise constant, on log haz scale
      lbound = negative_infinity();
    else
      lbound = 0;
    return lbound;
  }

  /**
  * Scale a vector of auxiliary parameters based on prior information
  *
  * @param aux_unscaled A vector, the unscaled auxiliary parameters
  * @param prior_dist Integer, the type of prior distribution
  * @param prior_mean,prior_scale Vectors, the mean and scale
  *   of the prior distribution
  * @return A vector, corresponding to the scaled auxiliary parameters
  */
  vector make_basehaz_coef(vector aux_unscaled, int prior_dist,
                           vector prior_mean, vector prior_scale) {
    vector[rows(aux_unscaled)] aux;
    if (prior_dist == 0) // none
      aux = aux_unscaled;
    else {
      aux = prior_scale .* aux_unscaled;
      if (prior_dist <= 2) // normal or student_t
        aux += prior_mean;
    }
    return aux;
  }

  /**
  * Log-prior for baseline hazard parameters
  *
  * @param aux_unscaled Vector (potentially of length 1) of unscaled
  *   auxiliary parameter(s)
  * @param dist Integer specifying the type of prior distribution
  * @param df Real specifying the df for the prior distribution
  * @return nothing
  */
  void basehaz_lp(vector aux_unscaled, int dist, vector df) {
    if (dist > 0) {
      if (dist == 1)
        target += normal_lpdf(aux_unscaled | 0, 1);
      else if (dist == 2)
        target += student_t_lpdf(aux_unscaled | df, 0, 1);
      else
        target += exponential_lpdf(aux_unscaled | 1);
    }
  }

  /**
  * Take the linear predictor and collapse across lower level
  * units of the grouping factor clustered within patients, using
  * the function specified by 'grp_assoc'
  *
  * @param eta The linear predictor evaluated for all the lower
  *   level units, having some length greater than N.
  * @param grp_idx An N-by-2 two dimensional array providing the
  *   beginning and ending index of the lower level units in eta that
  *   correspond to patient n (where n = 1,...,N).
  * @param grp_assoc The method for collapsing across the lower
  *   level units; 1=sum, 2=mean, 3=min, 4=max.
  * @return A vector
  */
  vector collapse_within_groups(vector eta, int[,] grp_idx,
                                int grp_assoc) {
    int N = size(grp_idx);
    vector[N] val;
    if (grp_assoc == 1) { // sum of lower level clusters
      for (n in 1:N)
        val[n] = sum(eta[grp_idx[n,1]:grp_idx[n,2]]);
    }
    else if (grp_assoc == 2) { // mean of lower level clusters
      for (n in 1:N)
        val[n] = mean(eta[grp_idx[n,1]:grp_idx[n,2]]);
    }
    else if (grp_assoc == 3) { // min of lower level clusters
      for (n in 1:N)
        val[n] = min(eta[grp_idx[n,1]:grp_idx[n,2]]);
    }
    else if (grp_assoc == 4) { // max of lower level clusters
      for (n in 1:N)
        val[n] = max(eta[grp_idx[n,1]:grp_idx[n,2]]);
    }
    return val;
  }

  /**
  * Create a design matrix for a shared random effects association
  * structure in the joint model
  *
  * @param b Vector of group-specific coefficients
  * @param l An integer array with the number of levels for the factor(s) on
  *   the RHS of each |
  * @param p An integer array with the number of variables on the LHS of each |
  * @param pmat A matrix with the number variables on the LHS of each | in each
  *   longitudinal submodel. The rows correspond to each |, meaning the separate
  *   equations for each grouping variable, and the columns correspond to each
  *   longitudinal submodel. If subject ID is the only grouping variable then the
  *   matrix will have one row. If the joint model only has one longitudinal
  *   submodel then the matrix will have one column.
  * @param Npat Integer specifying number of individuals represented
  *   in vector b
  * @param qnodes The number of quadrature nodes
  * @param which_b Integer array specifying the indices
  *   of the random effects to use in the association structure
  * @param sum_size_which_b Integer specifying total number of
  *   random effects that are to be used in the association structure
  * @param size_which_b Integer array specifying number of random effects from
  *   each long submodel that are to be used in the association structure
  * @param t_i Integer specifying the index of the grouping factor that
  *   corresponds to the patient-level
  * @param M An integer specifying the number of longitudinal submodels
  * @return A matrix with the desired random effects represented
  *   in columns, and the individuals on the rows; the matrix is
  *   repeated (qnodes + 1) times (bounded by rows)
  */
  matrix make_x_assoc_shared_b(
    vector b, int[] l, int[] p, int[,] pmat, int Npat, int qnodes,
    int[] which_b, int sum_size_which_b, int[] size_which_b, int t_i, int M) {
    int prior_shift; // num. ranefs prior to subject-specific ranefs
    int start_store;
    int end_store;
    matrix[Npat,sum_size_which_b] temp;
    matrix[(Npat*(qnodes+1)),sum_size_which_b] x_assoc_shared_b;
    if (t_i == 1) prior_shift = 0;
    else prior_shift = sum(l[1:(t_i-1)]);
    for (i in 1:Npat) {
      int mark;
      int start_collect;  // index start of subject-specific ranefs for patient
      mark = 1;
      start_collect = prior_shift + (i - 1) * p[t_i];
      for (m in 1:M) {
        if (size_which_b[m] > 0) {
          int shift;  // num. subject-specific ranefs in prior submodels
          int j_shift; // shift in indexing of which_b vector
          if (m == 1) {
            shift = 0;
            j_shift = 0;
          }
          else {
            shift = sum(pmat[t_i, 1:(m-1)]);
            j_shift = sum(size_which_b[1:(m-1)]);
          }
          for (j in 1:size_which_b[m]) {
            int item_collect;   // subject-specific ranefs to select for current submodel
            item_collect = start_collect + shift + which_b[(j_shift + j)];
            temp[i,mark] = b[item_collect];
            mark += 1;
          }
        }
      }
    }
    for (i in 1:(qnodes+1)) {
      start_store = (i - 1) * Npat + 1;
      end_store   = i * Npat;
      x_assoc_shared_b[start_store:end_store,] = temp;
    }
  return x_assoc_shared_b;
  }

  /**
  * Create a design matrix for a shared fixed + random effects association
  * structure in the joint model
  *
  * @param b Vector of group-specific coefficients
  * @param l An integer array with the number of levels for the factor(s) on
  *   the RHS of each |
  * @param p An integer array with the number of variables on the LHS of each |
  * @param pmat A matrix with the number variables on the LHS of each | in each
  *   longitudinal submodel. The rows correspond to each |, meaning the separate
  *   equations for each grouping variable, and the columns correspond to each
  *   longitudinal submodel. If subject ID is the only grouping variable then the
  *   matrix will have one row. If the joint model only has one longitudinal
  *   submodel then the matrix will have one column.
  * @param Npat Integer specifying number of individuals represented
  *   in vector b
  * @param qnodes The number of quadrature nodes
  * @param which_b Integer array specifying the indices
  *   of the random effects to use in the association structure
  * @param sum_size_which_b Integer specifying total number of
  *   random effects that are to be used in the association structure
  * @param size_which_b Integer array specifying number of random effects from
  *   each long submodel that are to be used in the association structure
  * @param t_i Integer specifying the index of the grouping factor that
  *   corresponds to the patient-level
  * @param M An integer specifying the number of longitudinal submodels
  * @return A matrix with the desired random effects represented
  *   in columns, and the individuals on the rows; the matrix is
  *   repeated (qnodes + 1) times (bounded by rows)
  */
  matrix make_x_assoc_shared_coef(
    vector b, vector beta, int[] KM, int M, int t_i,
    int[] l, int[] p, int[,] pmat, int Npat, int qnodes,
    int sum_size_which_coef, int[] size_which_coef,
    int[] which_coef_zindex, int[] which_coef_xindex,
    int[] has_intercept, int[] has_intercept_nob,
    int[] has_intercept_lob, int[] has_intercept_upb,
    real[] gamma_nob, real[] gamma_lob, real[] gamma_upb) {

    // in the loops below:
    //   t_i should only really ever equal 1 (since shared_coef association
    //       structure is not allowed if there is more than one clustering level)
    //   i = levels (ie, individuals)
    //   j = indices of the shared random effecs
    //   m = models

    int t_shift;  // skip over group-level coefficients for earlier grouping factors
    int start_store;
    int end_store;
    matrix[Npat,sum_size_which_coef] temp;
    matrix[(Npat*(qnodes+1)),sum_size_which_coef] x_assoc_shared_coef;
    if (t_i == 1) t_shift = 0;
    else t_shift = sum(l[1:(t_i-1)]);
    for (i in 1:Npat) {
      int mark;    // counter for looping over shared coefficients
      int i_shift; // skip over group-level coefficients for earlier levels
      mark = 1;
      i_shift = (i - 1) * p[t_i];
      for (m in 1:M) {
        if (size_which_coef[m] > 0) {  // if model has shared coefficients
          int j_shift;  // skip over elements of which_coef_zindex vector that are associated with earlier submodels
          int m_shift;  // skip over individual i's group-level coefficients for earlier submodels
          int shift_nb;
          int shift_lb;
          int shift_ub;
          int shift_beta;
          if (m == 1) {
            j_shift = 0; m_shift = 0; shift_nb = 0;
            shift_lb = 0; shift_ub = 0; shift_beta = 0;
          }
          else {
            j_shift = sum(size_which_coef[1:(m-1)]);
            m_shift = sum(pmat[t_i, 1:(m-1)]);
            shift_nb = sum(has_intercept_nob[1:(m-1)]);
            shift_lb = sum(has_intercept_lob[1:(m-1)]);
            shift_ub = sum(has_intercept_upb[1:(m-1)]);
            shift_beta = sum(KM[1:(m-1)]);
          }
          for (j in 1:size_which_coef[m]) {
            int b_collect;      // group-level coefficients to extract for current i, j, m
            int beta_collect_m; // within-submodel index of fixed effect coefficient to extract
            int beta_collect;   // overall index of fixed effect coefficient to extract
            real coef;
            b_collect = t_shift + i_shift + m_shift + which_coef_zindex[(j_shift + j)];
            beta_collect_m = which_coef_xindex[(j_shift + j)];
            beta_collect = shift_beta + beta_collect_m;
            coef = b[b_collect];  // start with group-level coefficient
            if ((has_intercept[m] == 1) && (beta_collect == 1)) {
              // collect intercept
              if (has_intercept_nob[m] == 1)
                coef += gamma_nob[sum(has_intercept_nob[1:m])];
              else if (has_intercept_lob[m] == 1)
                coef += gamma_lob[sum(has_intercept_lob[1:m])];
              else if (has_intercept_upb[m] == 1)
                coef += gamma_upb[sum(has_intercept_upb[1:m])];
            }
            else if (has_intercept[m] == 1) {
              // collect fixed effect whilst recognising intercept term
              // isn't in beta and correcting for that in the indexing
              coef += beta[(beta_collect - 1)];
            }
            else
              coef += beta[beta_collect];

            temp[i, mark] = coef;
            mark += 1;  // move to next shared coefficient for individual i
          }
        }
      }
    }

    // repeat the temp matrix qnodes times (ie, rbind)
    for (i in 1:(qnodes+1)) {
      start_store = (i - 1) * Npat + 1;
      end_store   = i * Npat;
      x_assoc_shared_coef[start_store:end_store, ] = temp;
    }
  return x_assoc_shared_coef;
  }

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
}

data {
  // declares: M, has_aux, has_weights, resp_type, intercept_type,
  //   yNobs, yNeta, yK, t, p, l, q, len_theta_L, bN1, bK1, bK1_len
  //   bK1_idx, bN2, bK2, bK2_len, bK2_idx
    // population level dimensions
  int<lower=1,upper=3> M; // num submodels with data (limit of 3)
  int<lower=0,upper=1> has_aux[3]; // has auxiliary param
  int<lower=0,upper=1> has_weights; // has observation weights
  int<lower=0,upper=2> resp_type[3]; // 1=real,2=integer,0=none
  int<lower=0,upper=3> intercept_type[3]; // 1=unbounded,2=lob,3=upb,0=none
  int<lower=0> yNobs[3]; // num observations
  int<lower=0> yNeta[3]; // required length of eta
  int<lower=0> yK[3]; // num predictors

  // group level dimensions, for decov prior
  int<lower=0> t;    // num. terms (maybe 0) with a | in the glmer formula
  int<lower=1> p[t]; // num. variables on the LHS of each |
  int<lower=1> l[t]; // num. levels for the factor(s) on the RHS of each |
  int<lower=0> q;    // conceptually equals \sum_{i=1}^t p_i \times l_i
  int<lower=0> len_theta_L; // length of the theta_L vector

  // group level dimensions, for lkj prior

    // group factor 1
    int<lower=0> bN1; // num groups
    int<lower=0> bK1; // total num params
    int<lower=0> bK1_len[3]; // num params in each submodel
    int<lower=0> bK1_idx[3,2]; // beg/end index for group params

    // group factor 2
    int<lower=0> bN2; // num groups
    int<lower=0> bK2; // total num params
    int<lower=0> bK2_len[3]; // num params in each submodel
    int<lower=0> bK2_idx[3,2]; // beg/end index for group params


      // population level data
  int<lower=0> yInt1[resp_type[1] == 2 ? yNobs[1] : 0]; // integer responses
  int<lower=0> yInt2[resp_type[2] == 2 ? yNobs[2] : 0];
  int<lower=0> yInt3[resp_type[3] == 2 ? yNobs[3] : 0];
  vector[resp_type[1] == 1 ? yNobs[1] : 0] yReal1; // real responses
  vector[resp_type[2] == 1 ? yNobs[2] : 0] yReal2;
  vector[resp_type[3] == 1 ? yNobs[3] : 0] yReal3;
  matrix[yNeta[1],yK[1]] yX1; // fe design matrix
  matrix[yNeta[2],yK[2]] yX2;
  matrix[yNeta[3],yK[3]] yX3;
  vector[yK[1]] yXbar1; // predictor means
  vector[yK[2]] yXbar2;
  vector[yK[3]] yXbar3;

  // family and link (determined by 'append_mvmer_famlink' R function)
  // 1 = gaussian
  // 2 = gamma
  // 3 = inverse gaussian
  // 4 = bernoulli
  // 5 = binomial (n>1)
  // 6 = poisson
  // 7 = negative binomial
  int<lower=0> family[M];
  int<lower=0> link[M]; // varies by family

  // group level data, group factor 1
  vector[bK1_len[1] > 0 ? yNeta[1] : 0] y1_Z1[bK1_len[1]]; // re design matrix
  vector[bK1_len[2] > 0 ? yNeta[2] : 0] y2_Z1[bK1_len[2]];
  vector[bK1_len[3] > 0 ? yNeta[3] : 0] y3_Z1[bK1_len[3]];
  int<lower=0> y1_Z1_id[bK1_len[1] > 0 ? yNeta[1] : 0]; // group indexing for y1_Z1
  int<lower=0> y2_Z1_id[bK1_len[2] > 0 ? yNeta[2] : 0]; // group indexing for y2_Z1
  int<lower=0> y3_Z1_id[bK1_len[3] > 0 ? yNeta[3] : 0]; // group indexing for y3_Z1

  // group level data, group factor 2
  vector[bK2_len[1] > 0 ? yNeta[1] : 0] y1_Z2[bK2_len[1]]; // re design matrix
  vector[bK2_len[2] > 0 ? yNeta[2] : 0] y2_Z2[bK2_len[2]];
  vector[bK2_len[3] > 0 ? yNeta[3] : 0] y3_Z2[bK2_len[3]];
  int<lower=0> y1_Z2_id[bK2_len[1] > 0 ? yNeta[1] : 0]; // group indexing for y1_Z2
  int<lower=0> y2_Z2_id[bK2_len[2] > 0 ? yNeta[2] : 0]; // group indexing for y2_Z2
  int<lower=0> y3_Z2_id[bK2_len[3] > 0 ? yNeta[3] : 0]; // group indexing for y3_Z2

  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus,
  //   5 = laplace, 6 = lasso, 7 = product_normal
  int<lower=0,upper=7> y_prior_dist[3];
  int<lower=0,upper=2> y_prior_dist_for_intercept[M];

  // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
  int<lower=0,upper=3> y_prior_dist_for_aux[M];

  // prior family: 1 = decov, 2 = lkj
  int<lower=1,upper=2> prior_dist_for_cov;

  // flag indicating whether to draw from the prior
  int<lower=0,upper=1> prior_PD;  // 1 = yes

  // declares:
  //   e_prior_dist{_for_intercept,_for_aux}
  //   qnodes
  //   len_{epts,qpts,ipts}
  //   epts,qpts,ipts
  //   qwts,iwts
  //   basehaz_type

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

  // declares:
  //   a_{K,xbar}
  //   a_prior_dist, assoc, assoc_uses, has_assoc
  //   {sum_}size_which_b, which_b_zindex
  //   {sum_}size_which_coef, which_coef_{zindex,xindex}
  //   {sum_,sum_size_}which_interactions
  //   y_nrow{_auc}_cpts, auc_{qnodes,qwts}
  //   a_K_data, has_grp, grp_assoc, idx_grp
  //   idx_cpts
  //   y{1,2,3}_x_        {eta,eps,auc}_cpts
  //   y{1,2,3}_z{1,2}_   {eta,eps,auc}_cpts
  //   y{1,2,3}_z{1,2}_id_{eta,eps,auc}_cpts



}
