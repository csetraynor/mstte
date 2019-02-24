
# Return design matrices for evaluating longitudinal submodel quantities
# at specified quadrature points/times
#
# @param data A data frame, the data for the longitudinal submodel.
# @param assoc A list with information about the association structure for
#   the one longitudinal submodel.
# @param y_mod A named list returned by a call to handle_y_mod (the
#   fit for a single longitudinal submodel)
# @param grp_stuff A list with information about any lower level grouping
#   factors that are clustered within patients and how to handle them in
#   the association structure.
# @param ids,times The subject IDs and times vectors that correspond to the
#   event and quadrature times at which the design matrices will
#   need to be evaluated for the association structure.
# @param id_var The name on the ID variable.
# @param time_var The name of the time variable.
# @param epsilon The half-width of the central difference used for
#   numerically calculating the derivative of the design matrix for slope
#   based association structures.
# @param auc_qnodes Integer specifying the number of GK quadrature nodes to
#   use in the integral/AUC based association structures.
# @return The list returned by make_assoc_parts.
handle_assocmod2 <- function(data, assoc_ids, assoc, y_mod, e_mod, grp_stuff, meta) {
  
  if (!requireNamespace("dplyr"))
    stop("the 'dplyr' package must be installed to use this function")
  
  if (!requireNamespace("data.table"))
    stop2("the 'data.table' package must be installed to use this function.")
  
  
  id_var   <- meta$id_var
  time_var <- meta$time_var
  
  assoc_obs <- data[[id_var]] %in% assoc_ids
  data <- data[assoc_obs, ]

  # before turning data into a data.table (for a rolling merge against
  # the quadrature points) we want to make sure that the data does not
  # include any NAs for the predictors or assoc formula variables
  tt <- attr(y_mod$terms, "term.labels")
  tt <- c(tt, uapply(assoc[["which_formulas"]], all.vars))
  fm <- reformulate(tt, response = NULL)
  df <- get_all_vars(fm, data)
  df <- df[complete.cases(df), , drop = FALSE]
  df[[id_var]] <- as.integer(as.factor(df[[id_var]])) # ensures smoothness in rolling merge
  
  # declare df as a data.table for merging with quadrature points
  dt <- prepare_data_table(df,
                           id_var   = id_var,
                           time_var = time_var,
                           grp_var  = grp_stuff$grp_var) # grp_var may be NULL

  # design matrices for calculating association structure based on
  # (possibly lagged) eta, slope, auc and any interactions with data
  args <- list(use_function = make_assoc_parts_for_stan,
               newdata      = dt,
               y_mod        = y_mod,
               grp_stuff    = grp_stuff,
               meta         = meta,
               assoc        = assoc,
               ids          = e_mod$cids_a,
               times        = e_mod$cpts_a)
  
  do.call(make_assoc_parts2, args)
}

# Function to construct quantities, primarily design matrices (x, Zt), that
# will be used to evaluate the longitudinal submodel contributions to the
# association structure in the event submodel. For example, the design matrices
# evaluated at the quadpoints, quadpoints + eps, lagged quadpoints, auc quadpoints,
# and so on. Exactly what quantities are returned depends on what is specified
# in the use_function argument.
#
# @param use_function The function to call which will return the design
#   matrices for eta, eps, lag, auc, etc. Generally either
#   'make_assoc_parts_for_stan' or 'pp_data'.
# @param newdata A model frame used for constructing the design matrices
# @param assoc A list with information about the association structure for
#   the one longitudinal submodel.
# @param grp_stuff A list with information about any lower level grouping
#   factors that are clustered within patients and how to handle them in
#   the association structure.
# @param ids,times The subject IDs and times vectors that correspond to the
#   event/censoring and quadrature times at which the design matrices will
#   need to be evaluated for the association structure.
# @param id_var The name on the ID variable.
# @param time_var The name of the time variable.
# @param epsilon The half-width of the central difference used for
#   numerically calculating the derivative of the design matrix for slope
#   based association structures.
# @param auc_qnodes Integer specifying the number of GK quadrature nodes to
#   use in the integral/AUC based association structures.
# @param ... Additional arguments passes to use_function
# @return A named list
make_assoc_parts2 <- function(use_function = make_assoc_parts_for_stan,
                              newdata, assoc, grp_stuff, meta_stuff,
                              ids, times, ...) {
  
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  
  id_var     <- meta_stuff$id_var
  time_var   <- meta_stuff$time_var
  epsilon    <- meta_stuff$epsilon
  auc_qnodes <- meta_stuff$auc_qnodes
  
  eps_uses_derivative_of_x <- TRUE # experimental
  
  # Apply lag
  lag <- assoc[["which_lag"]]
  if (!lag == 0)
    times <- set_lag(times, lag)
  
  # Broadcast ids and times if there is lower level clustering
  if (grp_stuff$has_grp) {
    # grps corresponding to each id
    grps <- as.vector(unlist(grp_stuff$grp_list[as.character(ids)]))
    # freq by which to expand each ids and times element
    freq_seq <- grp_stuff$grp_freq[as.character(ids)]
    # rep each patient id and prediction time the required num of times
    ids   <- rep(ids,   freq_seq)
    times <- rep(times, freq_seq)
    # indices for collapsing across clusters within patients
    grp_idx <- get_idx_array(freq_seq)
  } else grps <- grp_idx <- NULL
  
  # Identify row in longitudinal data closest to event time or quadrature point
  #   NB if the quadrature point is earlier than the first observation time,
  #   then covariates values are carried back to avoid missing values.
  #   In any other case, the observed covariates values from the most recent
  #   observation time preceeding the quadrature point are carried forward to
  #   represent the covariate value(s) at the quadrature point. (To avoid
  #   missingness there is no limit on how far forwards or how far backwards
  #   covariate values can be carried). If no time varying covariates are
  #   present in the longitudinal submodel (other than the time variable)
  #   then nothing is carried forward or backward.

  dataQ <- rolling_merge(data = newdata, ids = ids, times = times, grps = grps)
  
  dataQ <- handle_quadratic(dataQ, time_var)
  dataQ <- handle_cubic(dataQ, time_var)
  
  mod_eta <- use_function(newdata = dataQ, ...)
  
  # If association structure is based on slope, then calculate design
  # matrices under a time shift of epsilon
  sel_slope <- grep("etaslope", names(assoc))
  if (any(unlist(assoc[sel_slope]))) {
    if (eps_uses_derivative_of_x) {
      # slope is evaluated by passing Stan the derivatives of the X and Z
      # design matrices directly, each evaluated using central differences
      # with a half-width equal to epsilon
      dataQ_pos <- dataQ_neg <- dataQ
      dataQ_neg[[time_var]] <- dataQ_neg[[time_var]] - epsilon
      dataQ_pos[[time_var]] <- dataQ_pos[[time_var]] + epsilon
      mod_neg <- use_function(newdata = dataQ_neg, ...)
      mod_pos <- use_function(newdata = dataQ_pos, ...)
      mod_eps <- mod_pos
      mod_eps$x     <- (mod_pos$x     - mod_neg$x    ) / (2 * epsilon) # derivative of x
      mod_eps$xtemp <- (mod_pos$xtemp - mod_neg$xtemp) / (2 * epsilon) # derivative of centered x?
      mod_eps$z <- xapply(mod_pos$z, mod_neg$z,                        # derivative of z
                          FUN = function(x, y) (x - y) / (2 * epsilon))
      if (!is.null(mod_eps$Zt))
        mod_eps$Zt <- (mod_pos$Zt - mod_neg$Zt) / (2 * epsilon)
    } else {
      # slope is evaluated by passing Stan the X and Z design matrices under
      # a time shift of epsilon and then evaluating the derivative of the
      # linear predictor in Stan using a one-sided difference
      dataQ_eps <- dataQ
      dataQ_eps[[time_var]] <- dataQ_eps[[time_var]] + epsilon
      mod_eps <- use_function(newdata = dataQ_eps, ...)
    }
  } else mod_eps <- NULL
  
  # If association structure is based on area under the marker trajectory, then
  # calculate design matrices at the subquadrature points
  sel_auc <- grep("etaauc|muauc", names(assoc))
  if (any(unlist(assoc[sel_auc]))) {
    if (grp_stuff$has_grp)
      stop2("'etaauc' and 'muauc' not yet implemented when there is a grouping ",
            "factor clustered within patients.")
    # Return a design matrix that is (qnodes * auc_qnodes * Npat) rows
    auc_qpts <- uapply(times, function(x)
      lapply(get_quadpoints(auc_qnodes)$points, unstandardise_qpts, 0, x))
    auc_qwts <- uapply(times, function(x)
      lapply(get_quadpoints(auc_qnodes)$weights, unstandardise_qwts, 0, x))
    ids2 <- rep(ids, each = auc_qnodes)
    dataQ_auc <- rolling_merge(data = newdata, ids = ids2, times = auc_qpts)
    mod_auc <- use_function(newdata = dataQ_auc, ...)
  } else mod_auc <- auc_qpts <- auc_qwts <- NULL
  
  # If association structure is based on interactions with data, then calculate
  # the design matrix which will be multiplied by etavalue, etaslope, muvalue or muslope
  sel_data <- grep("_data", names(assoc), value = TRUE)
  X_data <- xapply(sel_data, FUN = function(i) {
    form <- assoc[["which_formulas"]][[i]]
    if (length(form)) {
      form <- as.formula(form)
      vars <- rownames(attr(terms.formula(form), "factors"))
      if (is.null(vars))
        stop2("No variables found in the formula for the '", i, "' association structure.")
      sel <- which(!vars %in% colnames(dataQ))
      if (length(sel))
        stop2("The following variables were specified in the formula for the '", i,
              "' association structure, but they cannot be found in the data: ",
              paste0(vars[sel], collapse = ", "))
      mf <- stats::model.frame(form, data = dataQ)
      X <- stats::model.matrix(form, data = mf)
      X <- drop_intercept(X)
      if (!ncol(X))
        stop2("Bug found: A formula was specified for the '", i, "' association ",
              "structure, but the resulting design matrix has no columns.")
    } else {
      X <- matrix(0, nrow(dataQ), 0)
    }
    X
  })
  K_data <- sapply(X_data, ncol)
  X_bind_data <- do.call(cbind, X_data)
  
  ret <- nlist(times, mod_eta, mod_eps, mod_auc, K_data, X_data, X_bind_data, grp_stuff)
  
  structure(ret,
            times      = times,
            lag        = lag,
            epsilon    = epsilon,
            grp_idx    = grp_idx,
            auc_qnodes = auc_qnodes,
            auc_qpts   = auc_qpts,
            auc_qwts   = auc_qwts,
            eps_uses_derivative_of_x = eps_uses_derivative_of_x)
}



# Check that the observation times for the longitudinal submodel are all
# positive and not observed after the individual's event time
#
# @param data A data frame (data for one longitudinal submodel)
# @param eventtimes A named numeric vector with the event time for each
#   individual. The vector names should be the individual ids.
# @param id_var,time_var The ID and time variable in the longitudinal data.
# @return Nothing.
validate_observation_times <-function(data, exittime, id_var, time_var) {
  if (!time_var %in% colnames(data))
    STOP_no_var(time_var)
  if (!id_var %in% colnames(data))
    STOP_no_var(id_var)
  if (any(data[[time_var]] < 0))
    stop2("Values for the time variable (", time_var, ") should not be negative.")
  mt  <- tapply(data[[time_var]], factor(data[[id_var]]), max) # max observation time
  nms <- names(mt)                                       # patient IDs
  if (is.null(nms))
    stop2("Bug found: cannot find names in the vector of exit times.")
  sel <- which(sapply(nms, FUN = function(i) mt[i] > exittime[i]))
  if (length(sel))
    stop2("The following individuals have observation times in the longitudinal ",
          "data that are later than their event time: ", comma(nms[sel]))
}

validate_observation_mstimes <- function(assoc_obs, mod, id_var, time_var, data){
  validate_observation_times(
    data = data[assoc_obs, ],
    exittime = mod$t_end_a,
    id_var     = id_var,
    time_var   = time_var
  )
}


# Validate the user input to the lag_assoc argument of stan_jm
#
# @param lag_assoc The user input to the lag_assoc argument
# @param M Integer specifying the number of longitudinal submodels
validate_lag_assoc <- function(lag_assoc, M) {
  if (length(lag_assoc) == 1L)
    lag_assoc <- rep(lag_assoc, M)
  if (!length(lag_assoc) == M)
    stop2("'lag_assoc' should length 1 or length equal to the ",
          "number of markers (", M, ").")
  if (!is.numeric(lag_assoc))
    stop2("'lag_assoc' must be numeric.")
  if (any(lag_assoc < 0))
    stop2("'lag_assoc' must be non-negative.")
  lag_assoc
}


# Filter associated times between the longitudinal model and hazard sumbmodel
#
# @param x Long dataset
# @param t_var time variable from long dataset.
# @param t time event vector
assoc_obstime <- function(e_mod, id_trans, long, meta){
  time_var = meta$time_var
  id_var = meta$id_var
  time_start = e_mod$t_start
  
  tmp <- data.frame(id = id_trans,
                    e_time = e_mod$t_end_a,
                    s_time = e_mod$t_start)
  colnames(tmp) <- c(id_var, "e_time", "s_time")
  
  tmpx <- dplyr::left_join(long, tmp, by = id_var)
  
  #out <- (tmpx[[time_var]] < tmpx[["e_time"]]) & (tmpx[[time_var]] >= tmpx[["s_time"]])
  
  out <- (tmpx[[time_var]] <= tmpx[["e_time"]]) 
  out[is.na(out)] <- FALSE
  return(out)
}


# Filter associated times between the longitudinal model and hazard sumbmodel
#
# @param x Long dataset
# @param t_var time variable from long dataset.
# @param t time event vector
match_obs <- function(newdataMs, id_trans, long, time_var, id_var, time_start_var){

  id_trans <- as.integer(id_trans)
  tmp <- dplyr::data_frame(id = id_trans,
                    e_time = newdataMs[[time_start_var]] + newdataMs[[time_var]],
                    s_time = newdataMs[[time_start_var]])
  colnames(tmp) <- c(id_var, "e_time", "s_time")
  
  tmpx <- dplyr::left_join(long, tmp, by = id_var)
  
  #out <- (tmpx[[time_var]] < tmpx[["e_time"]]) & (tmpx[[time_var]] >= tmpx[["s_time"]])
  
  out <- (tmpx[[time_var]] <= tmpx[["e_time"]]) 
  out[is.na(out)] <- FALSE
  return(out)
}