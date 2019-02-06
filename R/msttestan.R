# Function to create a stansurv object (fitted model object)
#
# @param object A list returned by a call to stan_surv
# @return A stansurv object
#
msttestan <- function(object) {

  alg     <- object$algorithm
  opt     <- alg == "optimizing"
  mcmc    <- alg == "sampling"
  stanfit <- object$stanfit
  basehaz <- object$basehaz
  K       <- lapply(object$x, function(X) NCOL(X) )

  if (opt)
    stop2("Optimisation not implemented for 'stanmstte' objects.")

  stan_summary <- make_stan_summary(stanfit)

  # number of parameters
  nvars  <- sum(
    nonrecuapply(seq_along(object$x), function(i)
      ncol(object$x[[i]]) + has_intercept(basehaz[[i]]) + basehaz[[i]]$nvars
    )
  )
  # obtain medians
  coefs     <- stan_summary[seq(nvars), select_median(alg)]
  coefs_nms <- rownames(stan_summary)[seq(nvars)]
  names(coefs) <- coefs_nms # ensure parameter names are retained

  # obtain standard errors and covariance matrix
  stanmat <- as.matrix(stanfit)[, seq(nvars), drop = FALSE]
  colnames(stanmat) <- coefs_nms
  ses <- apply(stanmat, 2L, mad)
  covmat <- cov(stanmat)

  # for mcmc only
  if (mcmc) {
    check_rhats(stan_summary[, "Rhat"])    # check rhats for all parameters
    runtime <- get_runtime(object$stanfit) # run time (in mins)
  }

  # NOT IMPLEMENTED QUADRATURE OR TDE.

  # return object of class 'stansurv'
  out <- nlist(
    coefficients  = coefs,
    ses           = ses,
    covmat        = covmat,
    formula       = object$formula,
    has_tde       = object$has_tde,
    has_quadrature= object$has_quadrature,
    terms         = object$terms,
    data          = object$data,
    transition_labels = object$transition_labels,
    model_frame   = object$model_frame,
    x             = object$x,
    s_cpts        = object$s_cpts,
    entrytime     = object$t_beg,
    eventtime     = object$t_end,
    event         = object$event,
    delayed       = object$delayed,
    basehaz       = object$basehaz,
    nobs          = object$nobs,
    nevents       = object$nevents,
    nlcens        = object$nlcens,
    nrcens        = object$nrcens,
    nicens        = object$nicens,
    ncensor       = object$ncensor,
    ndelayed      = object$ndelayed,
    qnodes        = object$qnodes,
    prior.info    = object$prior_info,
    algorithm     = object$algorithm,
    stan_function = object$stan_function,
    call          = object$call,
    runtime       = if (mcmc) runtime else NULL,
    rstan_version    = utils::packageVersion("rstan"),
    rstanarm_version = utils::packageVersion("rstanarm"),
    stan_summary,
    stanfit
  )
  out <- rm_null(out, recursive = FALSE)

  structure(out, class = c("stanmstte"))
}


#---------- internal

# Return the model fitting time in minutes.
#
# @param stanfit An object of class 'stanfit'.
# @return A matrix of runtimes, stratified by warmup/sampling and chain/overall.
get_runtime <- function(stanfit) {
  tt <- rstan::get_elapsed_time(stanfit)
  tt <- round(tt / 60, digits = 1L)    # time per chain
  tt <- cbind(tt, total = rowSums(tt)) # time per chain & overall
}
