# Construct a list with information about the multi-state submodel(s)
#
# @param formula The model formula for the event submodel
# @param data The data for the event submodel
# @param qnodes An integer specifying the number of GK quadrature nodes
# @param id_var The name of the ID variable
# @param y_id_list A character vector with a unique list of subject IDs
#   (factor levels) that appeared in the longitudinal submodels
# @return A named list with the following elements:
#   mod: The fitted Cox model.
#   entrytime: Named vector of numeric entry times.
#   eventtime: Named vector of numeric event times.
#   status: Named vector of event/failure indicators.
#   Npat: Number of individuals.
#   Nevents: Total number of events/failures.
#   id_list: A vector of unique subject IDs, as a factor.
#   qnodes: The number of GK quadrature nodes.
#   qwts,qpts: Vector of unstandardised quadrature weights and points.
#     The vector is ordered such that the first Npat items are the
#     weights/locations of the first quadrature point, then the second
#     Npat items are the weights/locations for the second quadrature
#     point, and so on.
#   qids: The subject IDs corresponding to each element of qwts/qpts.
#   epts: The event times, but only for individuals who were NOT censored
#     (i.e. those individual who had an event).
#   eids: The subject IDs corresponding to each element of epts.
#   cpts: Combined vector of failure and quadrature times: c(epts, qpts).
#   cids: Combined vector subject IDs: c(eids, qids).
#   Xq: The model matrix for the event submodel, centred and no intercept.
#   Xbar: Vector of column means for the event submodel model matrix.
#   K: Number of predictors for the event submodel.
#   norm_const: Scalar, the constant used to shift the event submodel
#     linear predictor (equal to the log of the mean incidence rate).
#   model_frame: The model frame for the fitted Cox model, but with the
#     subject ID variable also included.
#   tvc: Logical, if TRUE then a counting type Surv() object was used
#     in the fitted Cox model (ie. time varying covariates).
handle_ms_mod <- function(formula, data, basehaz, basehaz_ops = NULL, qnodes, t_start, meta) {
  
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function")
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")
  
  id_var      <- meta$id_var
  id_list     <- meta$id_list
  time_start_var  <- meta$time_start
  
  # parse formula, create model data & frame
  formula   <- parse_formula(formula, data)
  formula2  <- addto_formula(formula$formula, id_var) # includes id_var
  data      <- make_model_data (formula2, data)       # row subsetting etc.
  mf_stuff  <- make_model_frame(formula2, data)       # returns Surv object
  
  mf <- mf_stuff$mf # model frame
  mt <- mf_stuff$mt # model terms
  
  mf[[id_var]] <- promote_to_factor(mf[[id_var]]) # same as lme4
  ids <- factor(mf[[id_var]])
  
  # entry and exit times for each row of data
  t_beg <- make_t(mf, type = "beg") # entry time
  t_end <- make_t(mf, type = "end") # exit  time
  t_upp <- make_t(mf, type = "upp") # upper time for interval censoring
  
  # event indicator for each row of data
  status <- make_d(mf)
  
  event <- as.logical(status == 1)
  rcens <- as.logical(status == 0)
  lcens <- as.logical(status == 2)
  icens <- as.logical(status == 3)
  
  if (any(status < 0 || status > 3))
    stop2("Invalid status indicator in Surv object.")
  if (any(lcens))
    stop2("Cannot handle left censoring.")
  
  # delayed entry indicator for each row of data
  delayed  <- as.logical(!t_beg == 0)
  
  # time variables for stan
  t_event   <- t_end[event]   # exact event time
  t_rcens   <- t_end[rcens]   # right censoring time
  t_lcens   <- t_end[lcens]   # left  censoring time
  t_lower   <- t_end[icens]   # lower limit of interval censoring time
  t_upper   <- t_upp[icens]   # upper limit of interval censoring time
  t_delayed <- t_beg[delayed] # delayed entry time
  t_end_a   <- t_end + t_start
  # entry and exit times for each individual
  t_tmp <- t_end; t_tmp[icens] <- t_upp[icens]
  entrytime <- tapply(t_beg, ids, min)
  exittime  <- tapply(t_tmp, ids, max)
  
  # dimensions
  nevent   <- sum(event)
  nrcens   <- sum(rcens)
  nlcens   <- sum(lcens)
  nicens   <- sum(icens)
  ndelayed <- sum(delayed)
  
  # baseline hazard
  ok_basehaz <- c("weibull", "bs", "piecewise")
  ok_basehaz_ops <- get_ok_basehaz_ops(basehaz)
  basehaz <- handle_basehaz2(basehaz        = basehaz,
                             basehaz_ops    = basehaz_ops,
                             ok_basehaz     = ok_basehaz,
                             ok_basehaz_ops = ok_basehaz_ops,
                             times          = t_end,
                             status         = event,
                             upper_times    = t_upp)
  nvars <- basehaz$nvars # number of basehaz aux parameters
  
  # flag if intercept is required for baseline hazard
  has_intercept <- ai(has_intercept(basehaz))
  
  # standardised weights and nodes for quadrature
  qq <- get_quadpoints(nodes = qnodes)
  qp <- qq$points
  qw <- qq$weights
  
  # event times & ids (for events only)
  epts <- t_end[event] # event times
  epts_a <- t_end[event] + t_start[event]
  
  eids_a <- eids <- ids  [event] # subject ids
  
  # quadrature points & weights, evaluated for each row of data
  qpts <- uapply(qp, unstandardise_qpts, t_beg, t_end)
  qwts <- uapply(qw, unstandardise_qwts, t_beg, t_end)
  
  
  qpts_a <- uapply(qp, unstandardise_qpts, t_beg + t_start, t_end + t_start)
  qwts_a <- uapply(qw, unstandardise_qwts, t_beg + t_start, t_end + t_start)
  qids <- qids_a <- rep(ids, qnodes)
  
  # quadrature points & weights, evaluated at upper limit of rows w/ interval censoring
  if (nicens) {
    ipts <- uapply(qp, unstandardise_qpts, t_beg[icens], t_upper)
    iwts <- uapply(qw, unstandardise_qwts, t_beg[icens], t_upper)
    iids <- rep(ids[icens], qnodes)
    
    ipts_a <- uapply(qp, unstandardise_qpts, t_beg[icens] + t_start[icens], t_upper + t_start)
    iwts_a <- uapply(qw, unstandardise_qwts, t_beg[icens]+ t_start[icens], t_upper + t_start)
    iids_a <- rep(ids[icens], qnodes)
    
    
  } else {
    ipts <- rep(0,0)
    iwts <- rep(0,0)
    iids <- rep(0,0)
    ipts_a <- rep(0,0)
    iwts_a <- rep(0,0)
    iids_a <- rep(0,0)
  }
  
  cpts <- c(epts, qpts, ipts)
  cids <- c(eids, qids, iids)
  
  cpts_a <- c(epts_a, qpts_a, ipts_a)
  cids_a <- c(eids_a, qids_a, iids_a)
  
  # # quadrature points & weights, evaluated for rows with delayed entry
  # if (ndelayed) {
  #   qpts_delayed <- uapply(qp, unstandardise_qpts, 0, t_delayed) # qpts for entry time
  #   qwts_delayed <- uapply(qw, unstandardise_qwts, 0, t_delayed) # qwts for entry time
  # } else {
  #   qpts_delayed <- rep(0,0)
  #   qwts_delayed <- rep(0,0)
  # }
  
  # dimensions
  len_epts <- length(epts)
  len_qpts <- length(qpts)
  len_ipts <- length(ipts)
  len_cpts <- length(cpts)
  idx_cpts <- get_idx_array(c(len_epts, len_qpts, len_ipts))
  
  # basis terms for baseline hazard
  basis_epts <-  make_basis(epts, basehaz)
  basis_qpts <-  make_basis(qpts, basehaz)
  basis_ipts <-  make_basis(ipts, basehaz)
  
  # predictor matrices
  x <- make_x(formula$tf_form, mf)$x
  K <- ncol(x)
  
  x_cpts <- rbind(keep_rows(x, event),
                  rep_rows (x, times = qnodes),
                  rep_rows (keep_rows(x, icens), times = qnodes))
  
  # fit a cox model
  if (formula$surv_type %in% c("right", "counting")) {
    mod <- survival::coxph(formula$formula, data = data, x = TRUE)
  } else if (formula$surv_type %in% c("interval", "interval2")) {
    mod <- survival::survreg(formula$formula, data = data, x = TRUE)
  } else {
    stop("Bug found: Invalid Surv type.")
  }
  
  # calculate mean log incidence, used as a shift in log baseline hazard
  norm_const <- log(nevent / sum(t_end - t_beg))
  
  # get time start for initial, intermidiate and absorving states
  t_state = data[ ,time_start_var]
  
  nlist(mod,
        surv_type = formula$surv_type,
        qnodes,
        basehaz,
        has_intercept,
        has_icens = as.logical(nicens),
        model_frame = mf,
        entrytime,
        exittime,
        norm_const,
        t_beg,
        t_end,
        t_end_a,
        t_start,
        t_upp,
        t_event,
        t_rcens,
        t_lcens,
        t_lower,
        t_upper,
        t_state,
        status,
        nevent,
        nrcens,
        nlcens,
        nicens,
        ndelayed,
        len_epts,
        len_qpts,
        len_ipts,
        len_cpts,
        idx_cpts,
        epts,
        qpts,
        ipts,
        cpts,
        qwts,
        iwts,
        eids,
        qids,
        iids,
        cids,
        epts_a,
        qpts_a,
        ipts_a,
        cpts_a,
        qwts_a,
        iwts_a,
        eids_a,
        qids_a,
        iids_a,
        cids_a,
        basis_epts,
        basis_qpts,
        basis_ipts,
        x,
        x_cpts,
        x_bar = colMeans(x),
        K)
}


## Validate ids for ms joint model
## Check that in each transition there are the correct subjects.
# @param mod: ms_mod
# @param formula: msformula
# @param data: msdata
# @param id_ist: meta ids
validate_msjm_ids <- function(formula, data, id_list, meta){
  id_var <- meta$id_var
  # parse formula, create model data & frame
  formula   <- parse_formula(formula, data)
  formula2  <- addto_formula(formula$formula, id_var) # includes id_var
  data      <- make_model_data (formula2, data)       # row subsetting etc.
  mf_stuff  <- make_model_frame(formula2, data)       # returns Surv object
  
  mf <- mf_stuff$mf # model frame
  mt <- mf_stuff$mt # model terms
  
  mf[[id_var]] <- promote_to_factor(mf[[id_var]]) # same as lme4
  ids <- factor(mf[[id_var]])
  if(!assertthat::assert_that(all(id_list == ids)) ){
    stop2("The patient IDs (levels of the grouping factor) included in the longitudinal and event submodels do not match ")
  }
}


extract_id <- function(x, id_var){
  x[[id_var]]
}
