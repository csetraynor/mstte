
# Draw new group-specific parameters
#
# Run a Metropolis-Hastings algorithm to draw group-specific parameters for new
# groups conditional on new outcome data provided by the user. These parameters
# are required for the so-called "dynamic predictions" relevant to joint modelling
# of longitudinal and time-to-event data, whereby we wish to draw new group-specific
# parameters that condition on longitudinal data observed up to the current time t.
#
# @param object A stanjm object.
# @param stanmat Matrix of draws that are being used to generate the predictions.
# @param ndL A list of data frames with each element containing the prediction data 
#   for one longitudinal submodel.
# @param ndE A data frame with the prediction data for the event submodel.
# @param ids A vector of unique IDs for the individuals in the prediction data.
# @param times A vector of last known survival times for the individuals in the 
#   prediction data.
simulate_b_pars_ms <- function(object, stanmat, ndL, ndMS, ids, times, scale = 1.5) {
  
  # Preliminaries and dimensions
  p <- .p(object) # num of b pars for each grouping factor
  has_two_grp_factors <- (length(object$cnms) > 1L)  
  if (!has_two_grp_factors) { # one grouping factor
    b1_var <- object$id_var
    b1_p <- p[[b1_var]] # num of b pars for ID grouping factor
  } else { # more than one grouping factor
    if (get_M(object) > 1)
      STOP_dynpred("multivariate joint models with more than one grouping factor.")
    if (length(p) > 2L)
      STOP_dynpred("models with more than two grouping factors.")
    b1_var <- object$id_var
    b2_var <- grep(utils::glob2rx(b1_var), names(p), value = TRUE, invert = TRUE)
    b1_p <- p[[b1_var]] # num of b pars for ID grouping factor
    b2_p <- p[[b2_var]] # num of b pars for second grouping factor
    b2_n <- tapply(ndL[[1]][[b2_var]], 
                   ndL[[1]][[b1_var]], 
                   n_distinct) # num of unique levels for b2 within each ID
  }
  
  # Obtain a list with the posterior means for each parameter
  pars_means <- extract_pars(object, means = TRUE) 
  ids <- ids[[1]]
  ndE <- ndMS[[1]]
  times <- times[[1]]
  # Simulate new b pars
  cat("Drawing new random effects for", length(ids), "individuals. ")
  cat("Monitoring progress:\n")
  pb <- utils::txtProgressBar(min = 0, max = length(ids), style = 3)
  acceptance_rate <- c()
  b_new <- list()
  for (i in 1:length(ids)) {
    if (!has_two_grp_factors) { # one grouping factor
      len <- b1_p
    } else { # more than one grouping factor
      len <- b1_p + b2_p * b2_n[ids[[i]]] 
    }
    mat <- matrix(NA, nrow(stanmat), len)
    
    # Design matrices for individual i only
    dat_i <- .pp_data_jm(object, ndL, ndE, etimes = times[[i]], ids = ids[[i]])
    if (has_two_grp_factors) {
      dat_i$Ni <- b2_n[ids[[i]]]
    }
    # Obtain mode and var-cov matrix of posterior distribution of new b pars
    # based on asymptotic assumptions, used as center and width of proposal
    # distribution in MH algorithm
    inits <- rep(0, len)
    val <- optim(inits, optim_fn, object = object, data = dat_i, 
                 pars = pars_means, method = "BFGS", hessian = TRUE)
    mu_i <- val$par                       # asymptotic mode of posterior
    sigma_i <- scale * solve(val$hessian) # (scaled) asymptotic vcov of posterior
    
    # Run MH algorithm for each individual
    b_current <- mu_i # asympotic mode used as init value for MH algorithm
    accept <- c()
    for (s in 1:nrow(stanmat)) {
      pars_s <- extract_pars(object, stanmat[s, , drop = FALSE])
      b_step <- mh_step(b_old = b_current, mu = mu_i, sigma = sigma_i, 
                        df = 4, object = object, data = dat_i, pars = pars_s)
      accept[s] <- any(!b_step == b_current)
      mat[s,] <- b_current <- b_step
    }
    new_nms <- unlist(sapply(dat_i$assoc_parts, function(x) x$mod_eta$Z_names))
    colnames(mat) <- paste0("b[", new_nms, "]")
    utils::setTxtProgressBar(pb, i)
    acceptance_rate[[paste0(object$id_var, ":", ids[i])]] <- mean(accept)
    b_new[[i]] <- mat
  }
  close(pb)
  
  # return stanmat with only the new b pars included
  b_new <- do.call("cbind", b_new)     # cbind new b pars for all individuals
  sel <- b_names(colnames(stanmat))    # stanmat cols containing old b pars
  stanmat <- stanmat[, -sel, drop = F] # drop old b pars from stanmat
  stanmat <- cbind(stanmat, b_new)     # add new b pars to stanmat
  structure(stanmat, b_new = b_new, acceptance_rate = acceptance_rate)
}

# The function to optimise, in order to obtain the asymptotic mode and var-cov 
# matrix of the posterior distribution for the new b pars
# 
# @param b The vector of b parameters
# @param object A stanjm object
# @param data Output from .pp_data_jm
# @param pars Output from extract_pars
optim_fn <- function(b, object, data, pars) {
  nms <- lapply(data$assoc_parts, function(x) x$mod_eta$Z_names)
  pars <- substitute_b_pars(object, data, pars, new_b = b, new_Z_names = nms)
  ll <- .ll_jm(object, data, pars, include_b = TRUE)
  return(-ll) # optimise -ll for full joint model 
}    

# Perform one iteration of the Metropolis-Hastings algorithm
# 
# @param b_old The current vector of b parameters
# @param mu The mean vector for the proposal distribution
# @param sigma The variance-covariance matrix for the proposal distribution
# @param object A stanjm object
# @param data Output from .pp_data_jm
# @param pars Output from extract_pars
mh_step <- function(b_old, mu, sigma, df, object, data, pars) {
  # New proposal for b vector
  b_new <- rmt(mu = mu, Sigma = sigma, df = df)
  # Calculate density for proposal distribution
  propdens_old <- dmt(x = b_old, mu = mu, Sigma = sigma, df = df)
  propdens_new <- dmt(x = b_new, mu = mu, Sigma = sigma, df = df)
  # Calculate density for target distribution
  nms <- lapply(data$assoc_parts, function(x) x$mod_eta$Z_names)
  pars_old <- substitute_b_pars(object, data, pars, new_b = b_old, new_Z_names = nms)
  pars_new <- substitute_b_pars(object, data, pars, new_b = b_new, new_Z_names = nms)
  targdens_old <- .ll_jm(object, data, pars_old, include_b = TRUE)
  targdens_new <- .ll_jm(object, data, pars_new, include_b = TRUE)
  # MH accept/reject step
  accept_ratio <- exp(targdens_new - targdens_old - propdens_new + propdens_old)
  if (accept_ratio >= runif(1)) return(b_new) else return(b_old)
}

# Function to add new b parameters to the stanmat
#
# @param object A stanjm object
# @param data Output from .pp_data_jm
# @param pars Output from extract_pars
# @param new_b A vector of new b pars, or a list of vectors, with 
#   each element being the new b pars for a single submodel.
# @param new_Z_names A vector, or a list of vectors, with the names 
#   for the new b pars.
substitute_b_pars <- function(object, data, pars, new_b, new_Z_names) {
  M <- get_M(object)
  if (!is(new_b, "list")) { # split b into submodels
    if (M == 1) {
      new_b <- list(new_b)
    } else {
      y_cnms  <- fetch(object$glmod, "z", "group_cnms")
      len_b <- sapply(y_cnms, function(x) length(unlist(x)))
      new_b <- split(new_b, rep(1:length(len_b), len_b))
    }
  }
  if (!is(new_Z_names, "list")) { # split Z_names into submodels
    if (M == 1) {
      new_b <- list(new_b)
    } else {
      y_cnms  <- fetch(object$glmod, "z", "group_cnms")
      len_b <- sapply(y_cnms, function(x) length(unlist(x)))
      new_Z_names <- split(new_Z_names, rep(1:length(len_b), len_b))
    }
  }  
  mapply(function(x, y) {
    if (!identical(is.vector(x), is.vector(y)))
      stop("Bug found: new_b and new_Z_names should both be vectors or lists of vectors.")
    if (!identical(length(x), length(y)))
      stop("Bug found: new_b and new_Z_names should be the same length.") 
  }, new_b, new_Z_names) 
  pars$b <- mapply(function(b, nms) {
    names(b) <- paste0("b[", nms, "]")
    t(b)
  }, new_b, new_Z_names, SIMPLIFY = FALSE)
  pars$stanmat <- pars$stanmat[, -b_names(colnames(pars$stanmat)), drop = FALSE]
  pars$stanmat <- do.call("cbind", c(list(pars$stanmat), pars$b))
  return(pars)
}

.p <- function(object) {
  sapply(object$cnms, length)
}


# Return the design matrices required for evaluating the linear predictor or
# log-likelihood in post-estimation functions for a \code{stan_jm} model
#
# @param object A stanmvreg object
# @param newdataLong A data frame or list of data frames with the new 
#   covariate data for the longitudinal submodel
# @param newdataEvent A data frame with the new covariate data for the
#   event submodel
# @param ids An optional vector of subject IDs specifying which individuals
#   should be included in the returned design matrices.
# @param etimes An optional vector of times at which the event submodel
#   design matrices should be evaluated (also used to determine the 
#   quadrature times). If NULL then times are taken to be the eventimes in
#   the fitted object (if newdataEvent is NULL) or in newdataEvent.
# @param long_parts,event_parts A logical specifying whether to return the
#   design matrices for the longitudinal and/or event submodels.
# @return A named list (with components M, Npat, ndL, ndE, yX, tZt, 
#   yZnames, eXq, assoc_parts) 
.pp_data_jm <- function(object, newdataLong = NULL, newdataEvent = NULL, 
                        ids = NULL, etimes = NULL, long_parts = TRUE, 
                        event_parts = TRUE) {
  M <- get_M(object)
  id_var   <- object$id_var
  time_var <- object$time_var
  
  if (!is.null(newdataLong) || !is.null(newdataEvent))
    newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
  
  # prediction data for longitudinal submodels
  ndL <- if (is.null(newdataLong)) 
    get_model_data(object)[1:M] else newdatas[1:M]
  
  # prediction data for event submodel
  ndE <- if (is.null(newdataEvent)) 
    get_model_data(object)[["Event"]] else newdatas[["Event"]]   
  
  # possibly subset
  if (!is.null(ids)) {
    ndL <- subset_ids(object, ndL, ids)
    ndE <- subset_ids(object, ndE, ids)
  }
  id_list <- unique(ndE[[id_var]]) # unique subject id list
  
  # evaluate the last known survival time and status
  if (!is.null(newdataEvent) && is.null(etimes)) {
    # prediction data for the event submodel was provided but  
    # no event times were explicitly specified by the user, so
    # they must be evaluated using the data frame
    surv <- eval(formula(object, m = "Event")[[2L]], ndE)
    etimes  <- unclass(surv)[,"time"]
    estatus <- unclass(surv)[,"status"]    
  } else if (is.null(etimes)) {
    # if no prediction data was provided then event times are 
    # taken from the fitted model
    etimes  <- object$eventtime[as.character(id_list)]
    estatus <- object$status[as.character(id_list)]
  } else { 
    # otherwise, event times ('etimes') are only directly specified for dynamic   
    # predictions via posterior_survfit in which case the 'etimes' correspond 
    # to the last known survival time and therefore we assume everyone has survived
    # up to that point (ie, set estatus = 0 for all individuals), this is true 
    # even if there is an event indicated in the data supplied by the user.
    estatus <- rep(0, length(etimes))
  }
  res <- nlist(M, Npat = length(id_list), ndL, ndE)
  
  if (long_parts && event_parts)
    lapply(ndL, function(x) {
      if (!time_var %in% colnames(x)) 
        STOP_no_var(time_var)
      if (!id_var %in% colnames(x)) 
        STOP_no_var(id_var)
      if (any(x[[time_var]] < 0))
        stop2("Values for the time variable (", time_var, ") should not be negative.")
      mt <- tapply(x[[time_var]], factor(x[[id_var]]), max)
      if (any(mt > etimes))
        stop2("There appears to be observation times in the longitudinal data that ",
              "are later than the event time specified in the 'etimes' argument.")      
    }) 
  
  # response and design matrices for longitudinal submodels
  if (long_parts) {
    y <- lapply(1:M, function(m) eval(formula(object, m = m)[[2L]], ndL[[m]]))
    ydat <- lapply(1:M, function(m) pp_data(object, ndL[[m]], m = m))
    yX <- fetch(ydat, "x")
    yZt <- fetch(ydat, "Zt")
    yZ_names <- fetch(ydat, "Z_names")
    flist <- lapply(ndL, function(x) factor(x[[id_var]]))
    res <- c(res, nlist(y, yX, yZt, yZ_names, flist))
  }
  
  # design matrices for event submodel and association structure
  if (event_parts) {
    qnodes <- object$qnodes
    qq <- get_quadpoints(qnodes)
    qtimes <- uapply(qq$points,  unstandardise_qpts, 0, etimes)
    qwts   <- uapply(qq$weights, unstandardise_qwts, 0, etimes)
    starttime <- deparse(formula(object, m = "Event")[[2L]][[2L]])
    edat <- prepare_data_table(ndE, id_var, time_var = starttime)
    id_rep <- rep(id_list, qnodes + 1)
    times <- c(etimes, qtimes) # times used to design event submodel matrices
    edat <- rolling_merge(edat, ids = id_rep, times = times)
    eXq  <- .pp_data_mer_x(object, newdata = edat, m = "Event")
    assoc_parts <- lapply(1:M, function(m) {
      ymf <- ndL[[m]]
      grp_stuff <- object$grp_stuff[[m]]
      if (grp_stuff$has_grp) {
        grp_stuff <- get_extra_grp_info( # update grp_info with new data
          grp_stuff, flist = ymf, id_var = id_var, 
          grp_assoc = grp_stuff$grp_assoc)
      }
      ymf <- prepare_data_table(ymf, id_var = id_var, time_var = time_var,
                                grp_var = grp_stuff$grp_var)
      make_assoc_parts(
        ymf, assoc = object$assoc[,m], id_var = id_var, time_var = time_var, 
        ids = id_rep, times = times, grp_stuff = grp_stuff,
        use_function = pp_data, object = object, m = m)
    })
    assoc_attr <- nlist(.Data = assoc_parts, qnodes, qtimes, qwts, etimes, estatus)
    assoc_parts <- do.call("structure", assoc_attr)
    res <- c(res, nlist(eXq, assoc_parts))
  }
  
  return(res)
}

# the functions below are heavily based on a combination of 
# lme4:::predict.merMod and lme4:::mkNewReTrms, although they do also have 
# substantial modifications
.pp_data_mer_x <- function(object, newdata, m = NULL, ...) {
  x <- get_x(object, m = m)
  if (is.null(newdata)) return(x)
  form <- if (is.null(m)) attr(object$glmod$fr, "formula") else 
    formula(object, m = m)
  L <- length(form)
  form[[L]] <- lme4::nobars(form[[L]])
  RHS <- formula(substitute(~R, list(R = form[[L]])))
  Terms <- terms(object, m = m)
  mf <- model.frame(object, m = m)
  ff <- formula(form)
  vars <- rownames(attr(terms.formula(ff), "factors"))
  mf <- mf[vars]
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  mfnew <- model.frame(delete.response(Terms), newdata, xlev = orig_levs)
  x <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(x, "contrasts"))
  return(x)
}
# Draw new group-specific parameters
#
# Run a Metropolis-Hastings algorithm to draw group-specific parameters for new
# groups conditional on new outcome data provided by the user. These parameters
# are required for the so-called "dynamic predictions" relevant to joint modelling
# of longitudinal and time-to-event data, whereby we wish to draw new group-specific
# parameters that condition on longitudinal data observed up to the current time t.
#
# @param object A stanjm object.
# @param stanmat Matrix of draws that are being used to generate the predictions.
# @param ndL A list of data frames with each element containing the prediction data 
#   for one longitudinal submodel.
# @param ndE A data frame with the prediction data for the event submodel.
# @param ids A vector of unique IDs for the individuals in the prediction data.
# @param times A vector of last known survival times for the individuals in the 
#   prediction data.
simulate_b_pars <- function(object, stanmat, ndL, ndE, ids, times, scale = 1.5) {
  
  # Preliminaries and dimensions
  p <- .p(object) # num of b pars for each grouping factor
  has_two_grp_factors <- (length(object$cnms) > 1L)  
  if (!has_two_grp_factors) { # one grouping factor
    b1_var <- object$id_var
    b1_p <- p[[b1_var]] # num of b pars for ID grouping factor
  } else { # more than one grouping factor
    if (get_M(object) > 1)
      STOP_dynpred("multivariate joint models with more than one grouping factor.")
    if (length(p) > 2L)
      STOP_dynpred("models with more than two grouping factors.")
    b1_var <- object$id_var
    b2_var <- grep(utils::glob2rx(b1_var), names(p), value = TRUE, invert = TRUE)
    b1_p <- p[[b1_var]] # num of b pars for ID grouping factor
    b2_p <- p[[b2_var]] # num of b pars for second grouping factor
    b2_n <- tapply(ndL[[1]][[b2_var]], 
                   ndL[[1]][[b1_var]], 
                   n_distinct) # num of unique levels for b2 within each ID
  }
  
  # Obtain a list with the posterior means for each parameter
  pars_means <- extract_pars(object, means = TRUE) 
  
  # Simulate new b pars
  cat("Drawing new random effects for", length(ids), "individuals. ")
  cat("Monitoring progress:\n")
  pb <- utils::txtProgressBar(min = 0, max = length(ids), style = 3)
  acceptance_rate <- c()
  b_new <- list()
  for (i in 1:length(ids)) {
    if (!has_two_grp_factors) { # one grouping factor
      len <- b1_p
    } else { # more than one grouping factor
      len <- b1_p + b2_p * b2_n[ids[[i]]] 
    }
    mat <- matrix(NA, nrow(stanmat), len)
    # Design matrices for individual i only
    dat_i <- .pp_data_jm(object, ndL, ndE, etimes = times[[i]], ids = ids[[i]])
    if (has_two_grp_factors) {
      dat_i$Ni <- b2_n[ids[[i]]]
    }
    # Obtain mode and var-cov matrix of posterior distribution of new b pars
    # based on asymptotic assumptions, used as center and width of proposal
    # distribution in MH algorithm
    inits <- rep(0, len)
    val <- optim(inits, optim_fn, object = object, data = dat_i, 
                 pars = pars_means, method = "BFGS", hessian = TRUE)
    mu_i <- val$par                       # asymptotic mode of posterior
    sigma_i <- scale * solve(val$hessian) # (scaled) asymptotic vcov of posterior
    
    # Run MH algorithm for each individual
    b_current <- mu_i # asympotic mode used as init value for MH algorithm
    accept <- c()
    for (s in 1:nrow(stanmat)) {
      pars_s <- extract_pars(object, stanmat[s, , drop = FALSE])
      b_step <- mh_step(b_old = b_current, mu = mu_i, sigma = sigma_i, 
                        df = 4, object = object, data = dat_i, pars = pars_s)
      accept[s] <- any(!b_step == b_current)
      mat[s,] <- b_current <- b_step
    }
    new_nms <- unlist(sapply(dat_i$assoc_parts, function(x) x$mod_eta$Z_names))
    colnames(mat) <- paste0("b[", new_nms, "]")
    utils::setTxtProgressBar(pb, i)
    acceptance_rate[[paste0(object$id_var, ":", ids[i])]] <- mean(accept)
    b_new[[i]] <- mat
  }
  close(pb)
  
  # return stanmat with only the new b pars included
  b_new <- do.call("cbind", b_new)     # cbind new b pars for all individuals
  sel <- b_names(colnames(stanmat))    # stanmat cols containing old b pars
  stanmat <- stanmat[, -sel, drop = F] # drop old b pars from stanmat
  stanmat <- cbind(stanmat, b_new)     # add new b pars to stanmat
  structure(stanmat, b_new = b_new, acceptance_rate = acceptance_rate)
}
