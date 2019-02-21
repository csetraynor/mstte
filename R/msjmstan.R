# Internal function to create a msjmstan object (fitted model object)
#
# @param object A list returned by a call to any of: stan_msjm
# @return A stanmvreg object
#
msjmstan <- function(object) {
  
  opt        <- object$algorithm == "optimizing"
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
  nms <- c(stanfit@sim$fnames_oi)
  coefs <- list()
  ses <- list()
  
  # Coefs and SEs for longitudinal submodel(s)                    
  if (is_msjm) {
    y_coefs <- lapply(1:M, function(m)
      stan_summary[c(nms$y[[m]], nms$y_b[[m]]), select_median(object$algorithm)])
    y_stanmat <- lapply(1:M, function(m) 
      as.matrix(stanfit)[, c(nms$y[[m]], nms$y_b[[m]]), drop = FALSE])
    y_ses <- lapply(y_stanmat, function(m) apply(m, 2L, mad))
    y_covmat <- lapply(y_stanmat, cov)
    for (m in 1:M) {
      rownames(y_covmat[[m]]) <- colnames(y_covmat[[m]]) <- 
        rownames(stan_summary)[c(nms$y[[m]], nms$y_b[[m]])]
    }
    # Remove padding
    coefs[1:M] <- list_nms(lapply(y_coefs, unpad_reTrms.default), M, stub = stub)
    ses[1:M]   <- list_nms(lapply(y_ses, unpad_reTrms.default), M, stub = stub)
    
  # Coefs and SEs for event submodel    
    e_coefs <- stan_summary[c(nms$e, nms$a), select_median(object$algorithm)]        
    if (length(e_coefs) == 1L) 
      names(e_coefs) <- rownames(stan_summary)[c(nms$e, nms$a)[1L]]
    e_stanmat <- as.matrix(stanfit)[, c(nms$e, nms$a), drop = FALSE]
    e_ses <- apply(e_stanmat, 2L, mad)    
    e_covmat <- cov(e_stanmat)
    rownames(e_covmat) <- colnames(e_covmat) <- 
      rownames(stan_summary)[c(nms$e, nms$a)]
    coefs$Event <- e_coefs
    ses$Event <- e_ses
  }
  
  # Covariance matrix for fixed effects                    
  stanmat <- as.matrix(stanfit)[, c(nms$alpha, nms$beta), drop = FALSE]
  covmat <- cov(stanmat)
  
  if (object$algorithm == "sampling") { # for MCMC fits only
    # Check Rhats for all parameters
    check_rhats(stan_summary[, "Rhat"])    
    # Run time (mins)
    times <- round((rstan::get_elapsed_time(object$stanfit))/60, digits = 1)
    times <- cbind(times, total = rowSums(times))      
  } 
  
  out <- nlist(
    formula       = list_nms(object$formula, M, stub),
    terms         = list_nms(object$terms, M, stub),
    coefficients  = coefs, 
    ses           = ses,
    covmat        = covmat,
    prior.weights = object$weights, 
    prior.info    = object$prior.info,
    algorithm     = object$algorithm,
    call          = object$call,
    stan_function = object$stan_function,
    runtime       = if (object$algorithm == "sampling") times else NULL,
    rstan_version    = utils::packageVersion("rstan"),
    rstanarm_version = utils::packageVersion("rstanarm"),
    stan_summary, stanfit
  )
  if (is_mvmer) {
    out$cnms      <- object$cnms
    out$flevels   <- object$flevels
    out$n_markers <- object$M
    out$n_grps    <- object$n_grps
    out$n_yobs    <- list_nms(object$n_yobs, M, stub)
    out$family    <- list_nms(object$family, M, stub)
    out$glmod     <- list_nms(object$glmod, M, stub)
    out$data      <- if (!is_jm) list_nms(object$data, M, stub) else NULL
    classes <- c("stanmvreg", "stanreg", "lmerMod")
  }
  if (is_jm) {
    out$id_var    <- object$id_var
    out$time_var  <- object$time_var
    out$n_subjects<- object$n_subjects
    out$n_events  <- sum(object$survmod$status == 1)
    out$eventtime <- object$survmod$eventtime
    out$status    <- object$survmod$status > 0
    out$survmod   <- object$survmod
    out$qnodes    <- object$qnodes
    out$epsilon   <- object$epsilon    
    out$assoc     <- object$assoc
    out$assocmod  <- list_nms(object$assocmod, M, stub) 
    out$dataLong  <- list_nms(object$dataLong, M, stub) 
    out$dataEvent <- object$dataEvent
    out$grp_stuff <- object$grp_stuff
    out$fr        <- object$fr
    classes <- c("stanjm", "stanmvreg", "stanreg", "lmerMod")
  }
  out <- rm_null(out, recursive = FALSE)
  structure(out, class = classes)
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

