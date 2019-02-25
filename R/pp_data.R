pp_data <-
  function(object,
           newdata = NULL,
           re.form = NULL,
           offset = NULL,
           m = NULL,
           ...) {
    # validate_stanreg_object(object)

    .pp_data(object, newdata = newdata, offset = offset, ...)
  }

#-------------  for models without lme4 structure  -------------------------

.pp_data <- function(object, newdata = NULL, offset = NULL, ...) {
  if (is(object, "gamm4")) {
    requireNamespace("mgcv", quietly = TRUE)
    if (is.null(newdata))   x <- predict(object$jam, type = "lpmatrix")
    else x <- predict(object$jam, newdata = newdata, type = "lpmatrix")
    if (is.null(offset)) 
      offset <- object$offset %ORifNULL% rep(0, nrow(x))
    return(nlist(x, offset))
  }
  if (is.null(newdata)) {
    x <- get_x(object)
    if (is.null(offset)) {
      offset <- object$offset %ORifNULL% rep(0, nrow(x))
    }
    if (inherits(object, "betareg")) {
      return(nlist(x, offset, z_betareg = object$z))
    }
    
    return(nlist(x, offset))
  }
  
  offset <- .pp_data_offset(object, newdata, offset)
  Terms <- delete.response(terms(object))[[1]]
  
  m <- model.frame(Terms, newdata, xlev = object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses"))) 
    .checkMFClasses(cl, m)
  x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  if (is(object, "polr") && !is_scobit(object)) 
    x <- x[,colnames(x) != "(Intercept)", drop = FALSE]
  
  if (inherits(object, "betareg")) {
    mf <- model.frame(delete.response(object$terms$precision), 
                      data = newdata, na.action = object$na.action, 
                      xlev = object$levels$precision)
    z_betareg <- model.matrix(object$terms$precision, mf, contrasts = object$contrasts$precision)
    return(nlist(x, offset, z_betareg))
  }
  
  return(nlist(x, offset))
}


#------------------  for models fit using stan_mstte  -----------------------
# Return the design matrices required for evaluating the linear predictor or
# log-likelihood in post-estimation functions for a \code{stan_mstte} model

.pp_data_mstte <- function(object,
                           newdata = NULL,
                           times   = NULL,
                           at_quadpoints = FALSE,
                           ...) {
  
  formula <- object$formula
  basehaz <- object$basehaz
  
  if (is.null(newdata))
    newdata <- get_model_data(object)
  
  # flags
  has_tde        <- object$has_tde
  has_quadrature <- object$has_quadrature
  
  # define dimensions and times for quadrature
  if (any(has_quadrature) && at_quadpoints) {
    
    # if (is.null(times))
    #   stop("Bug found: 'times' must be specified.")
    #
    # # error check time variables
    # if (!length(times) == nrow(newdata))
    #   stop("Bug found: length of 'times' should equal number rows in the data.")
    #
    # # number of nodes
    # qnodes <- object$qnodes
    #
    # # standardised weights and nodes for quadrature
    # qq <- get_quadpoints(nodes = qnodes)
    # qp <- qq$points
    # qw <- qq$weights
    #
    # # quadrature points & weights, evaluated for each row of data
    # pts <- uapply(qp, unstandardise_qpts, 0, times)
    # wts <- uapply(qw, unstandardise_qwts, 0, times)
    #
    # # id vector for quadrature points
    # ids <- rep(seq_along(times), times = qnodes)
    #  not implemented
    
  } else { # predictions don't require quadrature
    
    pts    <- times
    wts    <- lapply(times, function(t) rep(NA, length(t)) )
    ids    <- if(!is.null(object$id_vars) ) object$id_vars else  lapply(object$model_frame, function(mf) {
      row.names(mf)
    })
  }
  
  # time-fixed predictor matrix
  tf_form <- lapply(formula, function(f) reformulate_rhs(rhs(f$tf_form)) )
  mf <- lapply(seq_len(object$n_trans), function(nt) make_model_frame(tf_form[[nt]], newdata[[nt]], check_constant = FALSE)$mf )
  x  <- lapply(seq_len(object$n_trans), function(nt){
    out <- make_x(tf_form[[nt]], mf[[nt]], xlevs= object$xlevs, check_constant = FALSE)$x
    if(ncol(out) > 0){
      colnames(out) <- append_trans(colnames(out), nt, object$transition_labels[nt])
    }
    out <- as.data.frame(out)
    tibble::rownames_to_column(out, var = "id_for_passing_to__")
  }  )
  x <-  suppressMessages( Reduce(full_join_NA, x) )
  
  if (any(has_quadrature) && at_quadpoints) {
    # x <- rep_rows(x, times = qnodes)
    # ni
  }
  
  # time-varying predictor matrix
  if (any(has_tde) ) { # model has tde not implmented
    # if (at_quadpoints) {
    #   # expand covariate data
    #   newdata <- rep_rows(newdata, times = qnodes)
    # }
    # if (all(is.na(pts))) {
    #   # temporary replacement to avoid error in creating spline basis
    #   pts_tmp <- rep(0, length(pts))
    # } else {
    #   # else use prediction times or quadrature points
    #   pts_tmp <- pts
    # }
    # s <- make_s(formula = object$formula$td_form,
    #             data    = newdata,
    #             times   = pts_tmp,
    #             xlevs   = object$xlevs)
    # if (all(is.na(pts))) {
    #   # if pts were all NA then replace the time-varying predictor
    #   # matrix with all NA, but retain appropriate dimensions
    #   s[] <- NaN
    # }
  } else { # model does not have tde
    s <- matrix(0, length(pts), 0)
  }
  
  # return object
  return(nlist(pts,
               wts,
               ids,
               x,
               s,
               has_quadrature,
               at_quadpoints,
               qnodes = object$qnodes))
}



#--------------------  for models fit using stan_jm  -----------------------

# Return the design matrices required for evaluating the linear predictor or
# log-likelihood in post-estimation functions for a \code{stan_jm} model.
# Code is twiked to fit various transitions.
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
.pp_data_jm_ms <- function(object, newdataLong = NULL, newdataEvent = NULL, 
                           ids = NULL, etimes = NULL, long_parts = TRUE, 
                           event_parts = TRUE) {
  M <- get_M(object)
  id_var   <- object$id_var
  time_var <- object$time_var
  
  # prediction data for longitudinal submodels
  ndL <- newdataLong
  
  # prediction data for event submodel
  ndE <- newdataEvent
  
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


#--------------------  for models fit using stan_jm  -----------------------

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
.pp_data_jm2 <- function(object, newdataLong = NULL, newdataEvent = NULL, 
                         ids = NULL, etimes = NULL, long_parts = TRUE, 
                         event_parts = TRUE) {
  M <- get_M(object)
  id_var   <- object$id_var
  time_var <- object$time_var
  
  # prediction data for longitudinal submodels
  ndL <- newdataLong
  
  # prediction data for event submodel
  ndE <- newdataEvent
  
  # possibly subset
  if (!is.null(ids)) {
    ndL <- subset_ids(object, ndL, ids)[[1]]
    ndE <- subset_ids(object, ndE, ids)[[1]]
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
  ndL <- list(ndL)
  
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
    qnodes <- 15L
    qq <- get_quadpoints(qnodes)
    qtimes <- uapply(qq$points,  unstandardise_qpts, 0, etimes)
    qwts   <- uapply(qq$weights, unstandardise_qwts, 0, etimes)
  #  starttime <- deparse(formula(object)[[2L]])[[2L]]
    ndE[[id_var]] <- as.integer(as.factor(ndE[[id_var]]))
    edat <- prepare_data_table(ndE, id_var, time_var = object$time_var)
    id_rep <- rep(id_list, qnodes + 1)
    times <- c(etimes, qtimes) # times used to design event submodel matrices
    edat <- rolling_merge(edat, ids = id_rep, times = times)
    
    eXq  <- .pp_data_mer_x(object, newdata = edat, m = "Event", ndE)
    assoc_parts <- lapply(1:M, function(m) {
      ymf <- ndL[[m]]
      grp_stuff <- object$grp_stuff[[m]]
      if (grp_stuff$has_grp) {
        grp_stuff <- get_extra_grp_info( # update grp_info with new data
          grp_stuff, flist = ymf, id_var = id_var, 
          grp_assoc = grp_stuff$grp_assoc)
      }
      ymf[[id_var]] <- as.integer(as.factor(ymf[[id_var]]))
      ymf <- prepare_data_table(ymf, id_var = id_var, time_var = time_var,
                                grp_var = grp_stuff$grp_var)
      make_assoc_parts3(
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

null_or_zero <- function(x) {
  isTRUE(is.null(x) || all(x == 0))
}

.pp_data_offset <- function(object, newdata = NULL, offset = NULL) {
  if (is.null(newdata)) {
    # get offset from model object (should be null if no offset)
    if (is.null(offset)) 
      offset <- object$offset %ORifNULL% model.offset(model.frame(object))
  } else {
    if (!is.null(offset))
      stopifnot(length(offset) == nrow(newdata))
    else {
      # if newdata specified but not offset then confirm that model wasn't fit
   
      # with an offset (warning, not error)
      if (!is.null(object$call$offset) || 
          !null_or_zero(object$offset) ) {
          # !null_or_zero(model.offset(model.frame(object)))) {
        warning(
          "'offset' argument is NULL but it looks like you estimated ", 
          "the model using an offset term.", 
          call. = FALSE
        )
      }
      offset <- rep(0, nrow(newdata))
    }
  }
  return(offset)
}


# the functions below are heavily based on a combination of 
# lme4:::predict.merMod and lme4:::mkNewReTrms, although they do also have 
# substantial modifications
.pp_data_mer_x <- function(object, newdata, m = NULL, ndE, ...) {
  x <- get_x(object, m = m)
  if (is.null(newdata)) return(x)
  form <- if (is.null(m)) attr(object$glmod$fr, "formula") else 
    formula(object)[[2]]
  L <- length(form)
  form[[L]] <- lme4::nobars(form[[L]])
  RHS <- formula(substitute(~R, list(R = form[[L]])))
  Terms <- terms(object, m = m)

  mf <- ndE
  ff <- formula(form)
  vars <- rownames(attr(terms.formula(ff), "factors"))[-1]
  mf <- mf[vars]
  isFac <- vapply(mf, is.factor, FUN.VALUE = TRUE)
  isFac[attr(Terms, "response")] <- FALSE
  orig_levs <- if (length(isFac) == 0) 
    NULL else lapply(mf[isFac], levels)
  mfnew <- model.frame(delete.response(Terms[[2]]), newdata, xlev = orig_levs)
  x <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(x, "contrasts"))
  return(x)
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


