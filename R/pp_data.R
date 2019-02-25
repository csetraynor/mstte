
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

