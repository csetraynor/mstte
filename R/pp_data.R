
#------------------  for models fit using stan_idm  -----------------------

.pp_data_mstte <- function(object,
                         newdata = NULL,
                         times   = NULL,
                         at_quadpoints = FALSE,
                         k,
                         ...) {

  formula <- object$formula[[k]]
  basehaz <- object$basehaz[[k]]

  if (is.null(newdata))
    newdata <- get_model_data(object)

  #NOT IMPLEMENTED
  # flags
  has_tde        <- object$has_tde
  has_quadrature <- object$has_quadrature
  #
  # # define dimensions and times for quadrature
  # if (has_quadrature && at_quadpoints) {
  #
  #   if (is.null(times))
  #     stop("Bug found: 'times' must be specified.")
  #
  #   # error check time variables
  #   if (!length(times) == nrow(newdata))
  #     stop("Bug found: length of 'times' should equal number rows in the data.")
  #
  #   # number of nodes
  #   qnodes <- object$qnodes
  #
  #   # standardised weights and nodes for quadrature
  #   qq <- get_quadpoints(nodes = qnodes)
  #   qp <- qq$points
  #   qw <- qq$weights
  #
  #   # quadrature points & weights, evaluated for each row of data
  #   pts <- uapply(qp, unstandardise_qpts, 0, times)
  #   wts <- uapply(qw, unstandardise_qwts, 0, times)
  #
  #   # id vector for quadrature points
  #   ids <- rep(seq_along(times), times = qnodes)
  #
  # } else { # predictions don't require quadrature

  pts    <- times
  wts    <- rep(NA, length(times))
  ids    <- seq_along(times)

  # }

  # time-fixed predictor matrix
  tf_form <-  reformulate_rhs(rhs(formula$tf_form))
  mf <- make_model_frame(tf_form, newdata, check_constant = FALSE)$mf
  x  <-  make_x(tf_form, mf, xlevs= object$xlevs, check_constant = FALSE)$x

  # if (has_quadrature && at_quadpoints) {
  #   x <- rep_rows(x, times = qnodes)
  # }

  # time-varying predictor matrix
  # if (has_tde) { # model has tde  NOT IMPLEMENTED
  #   if (at_quadpoints) {
  #     # expand covariate data
  #     newdata <- rep_rows(newdata, times = qnodes)
  #   }
  #   if (all(is.na(pts))) {
  #     # temporary replacement to avoid error in creating spline basis
  #     pts_tmp <- rep(0, length(pts))
  #   } else {
  #     # else use prediction times or quadrature points
  #     pts_tmp <- pts
  #   }
  #   s <- make_s(formula = object$formula$td_form,
  #               data    = newdata,
  #               times   = pts_tmp,
  #               xlevs   = object$xlevs)
  #   if (all(is.na(pts))) {
  #     # if pts were all NA then replace the time-varying predictor
  #     # matrix with all NA, but retain appropriate dimensions
  #     s[] <- NaN
  #   }
  # } else { # model does not have tde
  s <- matrix(0, length(pts), 0)
  # }

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


.pp_data_surv <- function(object, newdata, times, type = "surv", ...) {
  dots <- list(...)
  h <- dots$h


  formula <- object$formula[[h]]
  basehaz <- object$basehaz[[h]]

  if (is.null(newdata))
    newdata <- object$model_data # --> need to create get_model_data method?...

  #----- check for delayed entry

  # throw error if delayed entry is present in prediction data
  if (uses.start.stop(object, h)) {
    check <- try({
      mf_tmp <- make_model_frame(formula$lhs_form, newdata)$mf
      t_beg <- make_t(mf_tmp, type = "beg")
      if (any(!t_beg == 0))
        stop2("'posterior_survfit' cannot handle non-zero start times.")
    })
  }

  # return status indicator if being used for log likelihood evaluation
  if (type == "ll") {
    mf_tmp <- make_model_frame(formula$lhs_form, newdata)$mf
    status <- make_d(mf_tmp)
  }

  # otherwise assume one row per individual, with no delayed entry
  # (i.e. all t_beg = 0)
  mf <- make_model_frame(formula$tf_form, newdata)$mf

  #----- define dimensions and times for quadrature

  # flags
  has_tde        <- object$has_tde[h]
  has_quadrature <- object$has_quadrature[h]

  if (has_quadrature) { # model uses quadrature

    # number of nodes
    qnodes <- object$qnodes

    # standardised weights and nodes for quadrature
    qq <- get_quadpoints(nodes = qnodes)
    qp <- qq$points
    qw <- qq$weights

    # quadrature points & weights, evaluated for each row of data
    qpts <- uapply(qp, unstandardise_qpts, 0, times) # qpts for exit time
    qwts <- uapply(qw, unstandardise_qwts, 0, times) # qwts for exit time

  }

  #----- basis terms for baseline hazard

  if (!has_quadrature) {
    basis  <- make_basis(times, basehaz)
    ibasis <- make_basis(times, basehaz, integrate = TRUE)
  } else {
    basis      <- make_basis(times, basehaz)
    basis_qpts <- make_basis(qpts, basehaz)
  }

  #----- predictor matrices

  # time-fixed predictor matrix
  x <- make_x(formula$tf_form, mf, xlevs = object$xlevs)$x

  if (has_quadrature)
    x_qpts <- rep_rows(x, times = qnodes)

  # time-varying predictor matrix
  if (has_tde) {
    s      <- make_s(formula = object$formula$td_form,
                     data    = newdata,
                     times   = times,
                     xlevs   = object$xlevs)
    s_qpts <- make_s(formula = object$formula$td_form,
                     data    = rep_rows(newdata, times = qnodes),
                     times   = qpts,
                     xlevs   = object$xlevs)
  } else if (has_quadrature) { # model does not have tde
    s      <- matrix(0,length(times),0)
    s_qpts <- matrix(0,length(qpts) ,0)
  }

  if (has_quadrature) {
    return(nlist(times,
                 qpts,
                 qwts,
                 x,
                 x_qpts,
                 s,
                 s_qpts,
                 basis,
                 basis_qpts,
                 has_quadrature,
                 status = if (type == "ll") status else NULL))
  } else {
    return(nlist(times,
                 x,
                 basis,
                 ibasis,
                 has_quadrature,
                 status = if (type == "ll") status else NULL))
  }
}
# Checks for the type of survival response, and whether the prediction
# call therefore requires an id variable identifying individuals
#
# @param object A stansurv object.
# @return A logical.
uses.tde        <- function(object) { object$has_tde }
uses.start.stop <- function(object, k) {
  object$formula[[k]]$surv_type == "counting" }
requires.idvar  <- function(object) { uses.start.stop(object) || uses.tde(object) }
