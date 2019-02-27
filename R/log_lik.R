# Evaluate the log baseline hazard at the specified times
# given the vector or matrix of MCMC draws for the baseline
# hazard coeffients / parameters
#
# @param log_haz A vector containing the log hazard for each
#   individual, evaluated at each of the quadrature points. The
#   vector should be ordered such that the first N elements contain
#   the log_haz evaluated for each individual at quadrature point 1,
#   then the next N elements are the log_haz evaluated for each
#   individual at quadrature point 2, and so on.
# @param qnodes Integer specifying the number of quadrature nodes
#   at which the log hazard was evaluated for each individual.
# @param qwts A vector of unstandardised GK quadrature weights.
# @return A vector or matrix of log survival probabilities.
evaluate_log_survival <- function(log_haz, qnodes, qwts) {
  UseMethod("evaluate_log_survival")
}

evaluate_log_survival.default <- function(log_haz, qnodes, qwts) {
  # convert log hazard to hazard
  haz <- exp(log_haz)
  # apply GK quadrature weights
  weighted_haz <- qwts * haz
  # sum quadrature points for each individual to get cum_haz
  splitting_vec <- rep(1:qnodes, each = length(haz) / qnodes)
  cum_haz <- Reduce('+', split(weighted_haz, splitting_vec))
  # return: -cum_haz == log survival probability
  -cum_haz
}

evaluate_log_survival.matrix <- function(log_haz, qnodes, qwts) {
  # convert log hazard to hazard
  haz <- exp(log_haz)
  # apply GK quadrature weights
  weighted_haz <- sweep(haz, 2L, qwts, `*`)
  # sum quadrature points for each individual to get cum_haz
  cum_haz <- Reduce('+', array2list(weighted_haz, nsplits = qnodes))
  # return: -cum_haz == log survival probability
  -cum_haz
}

# Evaluate the log baseline hazard at the specified times
# given the vector or matrix of MCMC draws for the baseline
# hazard coeffients / parameters
#
# @param times A vector of times.
# @param basehaz A list with info about the baseline hazard.
# @param coefs A vector or matrix of parameter estimates (MCMC draws).
# @return A vector or matrix, depending on the input type of coefs.
evaluate_log_basehaz2 <- function(times, basehaz, coefs) {
  type <- basehaz$type_name
  if (type == "weibull") {
    X  <- log(times) # log times
    B1 <- log(coefs) # log shape
    B2 <- coefs - 1  # shape - 1
    log_basehaz <- as.vector(B1) + linear_predictor(B2,X)
  } else if (type == "bs") {
    X <- predict(basehaz$bs_basis, times) # b-spline basis
    B <- coefs                            # b-spline coefs
    log_basehaz <- linear_predictor(B,X)
  } else {
    stop2("Not yet implemented for basehaz = ", type)
  }
  log_basehaz
}


evaluate_log_haz <- function(times, basehaz, betas, aux, intercept = NULL, x) {
  eta  <- linear_predictor(betas, x)
  args <- nlist(times, basehaz, aux,  intercept)
  do.call(evaluate_log_basehaz, args) + eta
}


evaluate_basehaz <- function(times, basehaz, aux, intercept = NULL) {
  exp(evaluate_log_basehaz(times = times, basehaz = basehaz,
                           aux = aux, intercept = intercept))
}

#-------------

# Evaluate the log baseline survival at the specified times given the
# vector or matrix of MCMC draws for the baseline hazard parameters
#
# @param times A vector of times.
# @param basehaz A list with info about the baseline hazard.
# @param aux,intercept A vector or matrix of parameter estimates (MCMC draws).
# @return A vector or matrix, depending on the input type of aux.
evaluate_log_basesurv <- function(times, basehaz, aux, intercept = NULL) {
  switch(get_basehaz_name(basehaz),
         "exp"       = log_basesurv_exponential(times, log_scale = intercept),
         "weibull"   = log_basesurv_weibull (times, shape = aux, log_scale = intercept),
         "gompertz"  = log_basesurv_gompertz(times, scale = aux, log_shape = intercept),
         "ms"        = log_basesurv_ms(times, coefs = aux, basis = basehaz$basis),
         stop2("Bug found: unknown type of baseline hazard."))
}

log_basesurv_exponential <- function(x, log_scale) {
  -linear_predictor(exp(log_scale), x)
}
log_basesurv_weibull  <- function(x, shape, log_scale) {
  -exp(as.vector(log_scale) + linear_predictor(shape, log(x)))
}
log_basesurv_gompertz <- function(x, scale, log_shape) {
  -(as.vector(exp(log_shape) / scale)) * (exp(linear_predictor(scale, x)) - 1)
}
log_basesurv_ms <- function(x, coefs, basis) {
  -linear_predictor(coefs, basis_matrix(x, basis = basis, integrate = TRUE))
}

evaluate_log_surv <- function(times, basehaz, betas, aux, intercept = NULL, x) {
  eta  <- linear_predictor(betas, x)
  args <- nlist(times, basehaz, aux,  intercept)
  c( do.call(evaluate_log_basesurv, args) ) * exp(eta)

}
