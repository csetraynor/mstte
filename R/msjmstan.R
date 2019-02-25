# Internal function to create a msjmstan object (fitted model object)
#
# @param object A list returned by a call to any of: stan_msjm
# @return A stanmvreg object
#
msjmstan <- function(object) {
  
  alg        <- object$algorithm
  mcmc       <- alg == "sampling"
  opt        <- alg == "optimizing"
  stanfit    <- object$stanfit
  M          <- object$M
  n_trans    <- object$n_trans
  transition_labels <- object$transition_labels
  is_mstte    <- is.stanms(object)
  is_msjm      <- is.stanmsjm(object)
  stub       <- if (is_msjm) "Long" else "y"
  
  if (opt)
    stop("Optimisation not implemented for stanmsjm objects.")
  
  
  stan_summary <- make_stan_summary(stanfit)
  
  # number of parameters
  nvars  <- nrow(stan_summary)
  
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

  structure( nlist(
    formula       = list_nms(ulist(object$formula), M, stub, transition_labels),
    terms         = list_nms(object$terms, M, stub, transition_labels),
    coefficients  = coefs, 
    ses           = ses,
    covmat        = covmat,
    prior.weights = object$weights, 
    prior.info    = object$prior.info,
    algorithm     = object$algorithm,
    call          = object$call,
    stan_function = object$stan_function,
    cnms      = object$cnms,
    flevels   = object$flevels,
    n_markers = object$M,
    time_var  = object$time_var,
    id_var    = object$id_var,
    basehaz   = object$basehaz,
    assoc     = object$assoc,
    id_list   = object$id_list,
    obs_list  = object$obs_list,
    n_grps    = object$n_grps,
    grp_stuff = object$grp_stuff,
    n_yobs    = list_nms(object$n_yobs, M, stub),
    family    = list_nms(object$family, M, stub),
    glmod     = list_nms(object$glmod, M, stub),
    dataLong  = object$dataLong, 
    dataMs    = object$dataMs,
    n_trans   = object$n_trans,
    transition_labels = object$transition_labels,
    ms_mod    = object$msmod,
    a_mod     = object$assocmod,
    glmod     = object$glmod,
    runtime       = if (object$algorithm == "sampling") runtime else NULL,
    rstan_version    = utils::packageVersion("rstan"),
   # mstte_version = utils::packageVersion("mstte"),
    stan_summary, stanfit
  ), class = "stanmsjm")
}



# Drop the extra reTrms from a matrix x
#
# @param x A matrix or array (e.g. the posterior sample or matrix of summary
#   stats)
# @param columns Do the columns (TRUE) or rows (FALSE) correspond to the 
#   variables?
unpad_reTrms <- function(x, ...) UseMethod("unpad_reTrms")
unpad_reTrms.default <- function(x, ...) {
  if (is.matrix(x) || is.array(x))
    return(unpad_reTrms.array(x, ...))
  keep <- !grepl("_NEW_", names(x), fixed = TRUE)
  x[keep]
}

unpad_reTrms.array <- function(x, columns = TRUE, ...) {
  ndim <- length(dim(x))
  if (ndim > 3)
    stop("'x' should be a matrix or 3-D array")
  
  nms <- if (columns) 
    last_dimnames(x) else rownames(x)
  keep <- !grepl("_NEW_", nms, fixed = TRUE)
  if (length(dim(x)) == 2) {
    x_keep <- if (columns) 
      x[, keep, drop = FALSE] else x[keep, , drop = FALSE]
  } else {
    x_keep <- if (columns) 
      x[, , keep, drop = FALSE] else x[keep, , , drop = FALSE]
  }
  return(x_keep)
}

make_b_nms <- function(group, m = NULL, stub = "Long") {
  group_nms <- names(group$cnms)
  b_nms <- character()
  m_stub <- if (!is.null(m)) get_m_stub(m, stub = stub) else NULL
  for (i in seq_along(group$cnms)) {
    nm <- group_nms[i]
    nms_i <- paste(group$cnms[[i]], nm)
    levels(group$flist[[nm]]) <- gsub(" ", "_", levels(group$flist[[nm]]))
    if (length(nms_i) == 1) {
      b_nms <- c(b_nms, paste0(m_stub, nms_i, ":", levels(group$flist[[nm]])))
    } else {
      b_nms <- c(b_nms, c(t(sapply(paste0(m_stub, nms_i), paste0, ":", 
                                   levels(group$flist[[nm]])))))
    }
  }
  return(b_nms)  
}

