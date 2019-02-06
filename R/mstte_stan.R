#' Bayesian multi-state models via Stan
#'
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Bayesian inference for multi-state models (sometimes known as models for
#' competing risk data). Currently, the command fits standard parametric
#' (exponential, Weibull and Gompertz) and flexible parametric (cubic
#' spline-based) hazard models on the hazard scale, with covariates included
#' under assumptions of either proportional or non-proportional hazards.
#'
#' @export
#' @importFrom splines bs
#' @import splines2
#'
#' @param formula A list of Surv formula
#' @param data A data frame containing the variables specified in
#'   \code{formula}.
#' @param basehaz A list character string indicating which baseline hazard to use
#'   for the event submodel. Current options are:
#'   \itemize{
#'     \item \code{"ms"}: a flexible parametric model using cubic M-splines to
#'     model the baseline hazard. The default locations for the internal knots,
#'     as well as the basis terms for the splines, are calculated with respect
#'     to time. If the model does \emph{not} include any time-dependendent
#'     effects then a closed form solution is available for both the hazard
#'     and cumulative hazard and so this approach should be relatively fast.
#'     On the other hand, if the model does include time-dependent effects then
#'     quadrature is used to evaluate the cumulative hazard at each MCMC
#'     iteration and, therefore, estimation of the model will be slower.
#'     \item \code{"bs"}: a flexible parametric model using cubic B-splines to
#'     model the \emph{log} baseline hazard. The default locations for the
#'     internal knots, as well as the basis terms for the splines, are calculated
#'     with respect to time. A closed form solution for the cumulative hazard
#'     is \strong{not} available regardless of whether or not the model includes
#'     time-dependent effects; instead, quadrature is used to evaluate
#'     the cumulative hazard at each MCMC iteration. Therefore, if your model
#'     does not include any time-dependent effects, then estimation using the
#'     \code{"ms"} baseline hazard will be faster.
#'     \item \code{"exp"}: an exponential distribution for the event times.
#'     (i.e. a constant baseline hazard)
#'     \item \code{"weibull"}: a Weibull distribution for the event times.
#'     \item \code{"gompertz"}: a Gompertz distribution for the event times.
#'   }
#' @param basehaz_ops A named list specifying options related to the baseline
#'   hazard. Currently this can include: \cr
#'   \itemize{
#'     \item \code{df}: a positive integer specifying the degrees of freedom
#'     for the M-splines or B-splines. An intercept is included in the spline
#'     basis and included in the count of the degrees of freedom, such that
#'     two boundary knots and \code{df - 4} internal knots are used to generate
#'     the cubic spline basis. The default is \code{df = 6}; that is, two
#'     boundary knots and two internal knots.
#'     \item \code{knots}: An optional numeric vector specifying internal
#'     knot locations for the M-splines or B-splines. Note that \code{knots}
#'     cannot be specified if \code{df} is specified. If \code{knots} are
#'     \strong{not} specified, then \code{df - 4} internal knots are placed
#'     at equally spaced percentiles of the distribution of uncensored event
#'     times.
#'   }
#'   Note that for the M-splines and B-splines - in addition to any internal
#'   \code{knots} - a lower boundary knot is placed at the earliest entry time
#'   and an upper boundary knot is placed at the latest event or censoring time.
#'   These boundary knot locations are the default and cannot be changed by the
#'   user.
#' @param qnodes The number of nodes to use for the Gauss-Kronrod quadrature
#'   that is used to evaluate the cumulative hazard when \code{basehaz = "bs"}
#'   or when time-dependent effects (i.e. non-proportional hazards) are
#'   specified. Options are 15 (the default), 11 or 7.
#' @param prior_intercept The prior distribution for the intercept. Note
#'   that there will only be an intercept parameter when \code{basehaz} is set
#'   equal to one of the standard parametric distributions, i.e. \code{"exp"},
#'   \code{"weibull"} or \code{"gompertz"}, in which case the intercept
#'   corresponds to the parameter \emph{log(lambda)} as defined in the
#'   \emph{stan_surv: Survival (Time-to-Event) Models} vignette. For the cubic
#'   spline-based baseline hazards there is no intercept parameter since it is
#'   absorbed into the spline basis and, therefore, the prior for the intercept
#'   is effectively specified as part of \code{prior_aux}.
#'
#'   Where relevant, \code{prior_intercept} can be a call to \code{normal},
#'   \code{student_t} or \code{cauchy}. See the \link[=priors]{priors help page}
#'   for details on these functions. Note however that default scale for
#'   \code{prior_intercept} is 20 for \code{stan_surv} models (rather than 10,
#'   which is the default scale used for \code{prior_intercept} by most
#'   \pkg{rstanarm} modelling functions). To omit a prior on the intercept
#'   ---i.e., to use a flat (improper) uniform prior--- \code{prior_intercept}
#'   can be set to \code{NULL}.
#' @param prior_aux The prior distribution for "auxiliary" parameters related to
#'   the baseline hazard. The relevant parameters differ depending
#'   on the type of baseline hazard specified in the \code{basehaz}
#'   argument. The following applies:
#'   \itemize{
#'     \item \code{basehaz = "ms"}: the auxiliary parameters are the coefficients
#'     for the M-spline basis terms on the baseline hazard. These parameters
#'     have a lower bound at zero.
#'     \item \code{basehaz = "bs"}: the auxiliary parameters are the coefficients
#'     for the B-spline basis terms on the log baseline hazard. These parameters
#'     are unbounded.
#'     \item \code{basehaz = "exp"}: there is \strong{no} auxiliary parameter,
#'     since the log scale parameter for the exponential distribution is
#'     incorporated as an intercept in the linear predictor.
#'     \item \code{basehaz = "weibull"}: the auxiliary parameter is the Weibull
#'     shape parameter, while the log scale parameter for the Weibull
#'     distribution is incorporated as an intercept in the linear predictor.
#'     The auxiliary parameter has a lower bound at zero.
#'     \item \code{basehaz = "gompertz"}: the auxiliary parameter is the Gompertz
#'     scale parameter, while the log shape parameter for the Gompertz
#'     distribution is incorporated as an intercept in the linear predictor.
#'     The auxiliary parameter has a lower bound at zero.
#'   }
#'   Currently, \code{prior_aux} can be a call to \code{normal}, \code{student_t}
#'   or \code{cauchy}. See \code{\link{priors}} for details on these functions.
#'   To omit a prior ---i.e., to use a flat (improper) uniform prior--- set
#'   \code{prior_aux} to \code{NULL}.
#' @param prior_smooth This is only relevant when time-dependent effects are
#'   specified in the model (i.e. the \code{tde()} function is used in the
#'   model formula. When that is the case, \code{prior_smooth} determines the
#'   prior distribution given to the hyperparameter (standard deviation)
#'   contained in a random-walk prior for the cubic B-spline coefficients used
#'   to model the time-dependent coefficient. Lower values for the hyperparameter
#'   yield a less a flexible smooth function for the time-dependent coefficient.
#'   Specifically, \code{prior_smooth} can be a call to \code{exponential} to
#'   use an exponential distribution, or \code{normal}, \code{student_t} or
#'   \code{cauchy}, which results in a half-normal, half-t, or half-Cauchy
#'   prior. See \code{\link{priors}} for details on these functions. To omit a
#'   prior ---i.e., to use a flat (improper) uniform prior--- set
#'   \code{prior_smooth} to \code{NULL}. The number of hyperparameters depends
#'   on the model specification (i.e. the number of time-dependent effects
#'   specified in the model) but a scalar prior will be recylced as necessary
#'   to the appropriate length.
#'
#' @details
#' \subsection{Time dependent effects (i.e. non-proportional hazards)}{
#'   By default, any covariate effects specified in the \code{formula} are
#'   included in the model under a proportional hazards assumption. To relax
#'   this assumption, it is possible to estimate a time-dependent coefficient
#'   for a given covariate. This can be specified in the model \code{formula}
#'   by wrapping the covariate name in the \code{tde()} function (note that
#'   this function is not an exported function, rather it is an internal function
#'   that can only be evaluated within the formula of a \code{stan_surv} call).
#'
#'   For example, if we wish to estimate a time-dependent effect for the
#'   covariate \code{sex} then we can specify \code{tde(sex)} in the
#'   \code{formula}, e.g. \code{Surv(time, status) ~ tde(sex) + age + trt}.
#'   The coefficient for \code{sex} will then be modelled
#'   using a flexible smooth function based on a cubic B-spline expansion of
#'   time.
#'
#'   The flexibility of the smooth function can be controlled in two ways:
#'   \itemize{
#'   \item First, through control of the prior distribution for the cubic B-spline
#'   coefficients that are used to model the time-dependent coefficient.
#'   Specifically, one can control the flexibility of the prior through
#'   the hyperparameter (standard deviation) of the random walk prior used
#'   for the B-spline coefficients; see the \code{prior_smooth} argument.
#'   \item Second, one can increase or decrease the number of degrees of
#'   freedom used for the cubic B-spline function that is used to model the
#'   time-dependent coefficient. By default the cubic B-spline basis is
#'   evaluated using 3 degrees of freedom (that is a cubic spline basis with
#'   boundary knots at the limits of the time range, but no internal knots).
#'   If you wish to increase the flexibility of the smooth function by using a
#'   greater number of degrees of freedom, then you can specify this as part
#'   of the \code{tde} function call in the model formula. For example, to
#'   use cubic B-splines with 7 degrees of freedom we could specify
#'   \code{tde(sex, df = 7)} in the model formula instead of just
#'   \code{tde(sex)}. See the \strong{Examples} section below for more
#'   details.
#'   }
#'   In practice, the default \code{tde()} function should provide sufficient
#'   flexibility for model most time-dependent effects. However, it is worth
#'   noting that the reliable estimation of a time-dependent effect usually
#'   requires a relatively large number of events in the data (e.g. >1000).
#' }
#'
#' @examples
#' \donttest{
#' #---------- Proportional hazards
#'
#' # Simulated data
#' library("SemiCompRisks")
#' library("scales")
#'
#'
#' # Data generation parameters
#' n <- 1500
#' beta1.true <- c(0.1, 0.1)
#' beta2.true <- c(0.2, 0.2)
#' beta3.true <- c(0.3, 0.3)
#' alpha1.true <- 0.12
#' alpha2.true <- 0.23
#' alpha3.true <- 0.34
#' kappa1.true <- 0.33
#' kappa2.true <- 0.11
#' kappa3.true <- 0.22
#' theta.true <- 0
#'
#' # Make design matrix with single binary covariate
#' x_c <- rbinom(n, size = 1, prob = 0.7)
#' x_m <- cbind(1, x_c)
#'
#' # Generate semicompeting data
#' dat_ID <- SemiCompRisks::simID(x1 = x_m, x2 = x_m, x3 = x_m,
#'                                beta1.true = beta1.true,
#'                                beta2.true = beta2.true,
#'                                beta3.true = beta3.true,
#'                                alpha1.true = alpha1.true,
#'                                alpha2.true = alpha2.true,
#'                                alpha3.true = alpha3.true,
#'                                kappa1.true = kappa1.true,
#'                                kappa2.true = kappa2.true,
#'                                kappa3.true = kappa3.true,
#'                                theta.true = theta.true,
#'                                cens = c(240, 360))
#'
#' dat_ID$x_c <- x_c
#' colnames(dat_ID)[1:4] <- c("R", "delta_R", "tt", "delta_T")
#'
#'
#'
#' formula01 <- Surv(time=rr,event=delta_R)~x_c
#' formula02 <- Surv(time=tt,event=delta_T)~x_c
#' formula12 <- Surv(time=difft,event=delta_T)~x_c
#'
#' }
mstte_stan <- function(formula = lapply(1:3, function (x)
  as.formula(Surv(time=time,event=status)~1) ),
                     transition_labels = NULL,
                     transmat,
                     prep = FALSE,
                     ids,
                     status,
                     times,
                     keep,
                     data,
                     basehaz = lapply(1:3, function(x)
                       "ms"),
                     basehaz_ops = NULL,
                     prior = lapply(1:3, function(x)
                       rstanarm::normal() ),
                     prior_intercept  = lapply(1:3, function(x)
                       rstanarm::normal() ),
                     prior_aux = lapply(1:3, function(x)
                       rstanarm::normal() ),
                     prior_PD        = FALSE,
                     algorithm       = c("sampling", "meanfield", "fullrank"),
                     adapt_delta     = 0.99, ...
){


  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------

  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function.")

  if(missing(formula)){
    stop2("A formula object has to be provided.")
  }

  nt <- length(formula) # number of transitions in the model

  if (missing(basehaz_ops))
    basehaz_ops <- NULL

  if (missing(data) || !inherits(data, "data.frame"))
    stop("'data' must be a data frame.")

  dots      <- list(...)
  algorithm <- match.arg(algorithm)

  ## Parse formula
  formula <- lapply(formula, function(f) parse_formula(f, data) )

  ## Create data
  data <-  lapply(unique(data$trans), function(i){
    data[data$trans == i,]
  })
  data <- lapply(seq_along(data), function(i) make_model_data(formula[[i]]$tf_form, data[[i]] ) )
  # row subsetting etc.

  pos <- sapply(data, function(d) nrow(d))

  #----------------
  # Construct data
  #----------------

  #----- model frame stuff

  mf_stuff <-  lapply(seq_along(formula), function(i) {
    make_model_frame(formula[[i]]$tf_form, data[[i]])
  })

  mf <- lapply(mf_stuff, function(m) m$mf)   # model frame
  mt <- lapply(mf_stuff, function(m) m$mt)  # model terms

  #----- dimensions and response vectors

  # entry and exit times for each row of data
  t_beg <- lapply(mf, function(m) make_t(m, type = "beg") ) # entry time
  t_end <- lapply(mf, function(m) make_t(m, type = "end") ) # exit time
  t_upp <- lapply(mf, function(m) make_t(m, type = "upp") ) # upper time for interval censoring

  # ensure no event or censoring times are zero (leads to degenerate
  # estimate for log hazard for most baseline hazards, due to log(0))

  for(i in seq_along(t_end ) ){
    check1 <- any(t_end[[i]] <= 0, na.rm = TRUE)
    if (check1)
      stop2("All event and censoring times must be greater than 0.")
  }

  # event indicator for each row of data
  status <- lapply(mf, function(m) make_d(m))

  for(i in seq_along(status)){
    if (any(status[[i]] < 0 || status[[i]] > 3))
      stop2("Invalid status indicator in formula.")
  }

  # delayed entry indicator for each row of data
  delayed <- lapply(t_beg, function(t)as.logical(!t == 0) )

  # time variables for stan
  t_event <- lapply(seq_along(t_end), function(t)
    t_end[[t]][status[[t]] == 1] )  # exact event time
  t_rcens <- lapply(seq_along(t_end), function(t)
    t_end[[t]][status[[t]] == 0] ) # right censoring time
  t_lcens <- lapply(seq_along(t_end),  function(t)
    t_end[[t]][status[[t]] == 2] ) # left censoring time
  t_icenl <- lapply(seq_along(t_end),  function(t)
    t_end[[t]][status[[t]] == 3] ) # lower limit of interval censoring time
  t_icenu <- lapply(seq_along(t_upp), function(t)
    t_upp[[t]][status[[t]] == 3] ) # upper limit of interval censoring time
  t_delay <- lapply(seq_along(t_beg),  function(t)
    t_beg[[t]][delayed[[t]]] )

  # calculate log crude event rate
  t_tmp <- lapply(seq_along(t_end), function(t) sum(rowMeans(cbind(t_end[[t]], t_upp[[t]] ), na.rm = TRUE) - t_beg[[t]] ))
  d_tmp <- lapply(status, function(s) sum(!s == 0) )
  log_crude_event_rate <- lapply(seq_along(t_tmp), function(t) log(d_tmp[[t]] / t_tmp[[t]]))

  # dimensions
  nevent <- lapply(status, function(s)  sum(s == 1))
  nrcens <- lapply(status,  function(s) sum(s == 0))
  nlcens <- lapply(status, function(s) sum(s == 2))
  nicens <- lapply(status, function(s) sum(s == 3))
  ndelay <- lapply(delayed, function(d) sum(d))

  #----- baseline hazard
  ok_basehaz <- c("exp", "weibull", "gompertz", "ms", "bs")
  ok_basehaz_ops <- lapply(basehaz, function(b) get_ok_basehaz_ops(b))

  basehaz <- lapply(seq_along(basehaz), function(b)
    SW( handle_basehaz_surv(basehaz = basehaz[[b]],
                            basehaz_ops = basehaz_ops[[b]],
                            ok_basehaz = ok_basehaz,
                            ok_basehaz_ops = ok_basehaz_ops[[b]],
                            times = t_end[[b]],
                            status = status[[b]],
                            min_t = min(t_beg[[b]]),
                            max_t = max(c(t_end[[b]], t_upp[[b]]), na.rm = TRUE ) )
    ))

  nvars <- lapply(basehaz, function(b) b$nvars)  # number of basehaz aux parameters
  # flag if intercept is required for baseline hazard
  has_intercept   <- lapply(basehaz, function(b) ai(has_intercept(b)) )
  has_quadrature <- lapply(basehaz,function(b) FALSE) #NOT IMPLEMENTED IN THIS RELEASE


  #----- basis terms for baseline hazard
  #       # NOT IMPLEMENTED
  #       #basis_cpts[[b]][[i]] <- make_basis(cpts[[b]], basehaz[[i]])


  basis_event <- lapply(seq_along(basehaz), function(i) make_basis(times = t_event[[i]], basehaz = basehaz[[i]]))
  ibasis_event <- lapply(seq_along(basehaz), function(i) make_basis(times = t_event[[i]], basehaz = basehaz[[i]], integrate = TRUE) )
  ibasis_rcens <- lapply(seq_along(basehaz), function(i) make_basis(times = t_rcens[[i]], basehaz = basehaz[[i]], integrate = TRUE) )
  ibasis_lcens <- lapply(seq_along(basehaz), function(i) make_basis(times = t_lcens[[i]], basehaz = basehaz[[i]], integrate = TRUE) )
  ibasis_icenl <- lapply(seq_along(basehaz), function(i) make_basis(times = t_icenl[[i]], basehaz = basehaz[[i]], integrate = TRUE) )
  ibasis_icenu <- lapply(seq_along(basehaz), function(i) make_basis(times = t_icenu[[i]], basehaz = basehaz[[i]], integrate = TRUE) )
  ibasis_delay <- lapply(seq_along(basehaz), function(i) make_basis(times = t_delay[[i]],basehaz =  basehaz[[i]], integrate = TRUE) )

  #----- predictor matrices

  # time-fixed predictor matrix
  x_stuff <- lapply(seq_along(mf), function(i) make_x(formula = formula[[i]]$tf_form, mf[[i]]))
  x   <- lapply(x_stuff,function(n) n$x)
  x_bar <- lapply(x_stuff, function(n) n$x_bar)
  # column means of predictor matrix

  x_centered <-  lapply(x_stuff, function(n) n$x_centered  )
  x_event <- lapply(seq_along(status), function(i) keep_rows(x_centered[[i]], status[[i]] == 1) )
  x_lcens <- lapply(seq_along(status), function(i) keep_rows(x_centered[[i]], status[[i]] == 2) )
  x_rcens <- lapply(seq_along(status), function(i) keep_rows(x_centered[[i]], status[[i]] == 0) )
  x_icens <- lapply(seq_along(status), function(i) keep_rows(x_centered[[i]], status[[i]] == 3) )
  x_delay <- lapply(seq_along(status),  function(i)keep_rows(x_centered[[i]], delayed[[i]]) )
  K <- lapply(x, function(i) ncol(i) )

  ## prepare data
  nevent = aau(nevent)
  nrcens = aau(nrcens)
  nlcens = aau(nlcens)
  nicens = aau(nicens)
  ndelay = aau(ndelay)

  posevent = cumsum(nevent)
  posrcens = cumsum(nrcens)
  poslcens = cumsum(nlcens)
  posnices = cumsum(nicens)
  posdelay = cumsum(ndelay)

  posbasis_event = aau(lapply(basis_event, function(b) nrow(b)))
  posibasis_event = aau(lapply(ibasis_event, function(b) nrow(b)))
  posibasis_rcens = aau(lapply(ibasis_rcens, function(b) nrow(b)))

  standata <- nlist(
    nt    = nt,
    s_K     = aau(K),
    nK    = sum(aau(K)),
    s_vars = aau(nvars),
    Nvars = sum(aau(nvars)),

    x_bar = aau(x_bar),

    has_intercept = aau(has_intercept),
    N_has_intercept = sum(aau(has_intercept)),

    type = aau(lapply(basehaz, function(b) b$type) ),

    log_crude_event_rate = aau(log_crude_event_rate),

    Nevent = sum(nevent),
    Nrcens = sum(nrcens),

    s_event = nevent,
    s_rcens = nrcens,

    t_event = aau(t_event),
    t_rcens = aau(t_rcens),

    x_event = aau(x_event),
    x_rcens = aau(x_rcens),

    Nxevent = length(aau(x_event)),
    Nxrcens = length(aau(x_rcens)),

    basis_event  = aau(basis_event),
    ibasis_event = aau(ibasis_event),
    ibasis_rcens = aau(ibasis_rcens),

    Nbasis_event  = length(aau(basis_event)),
    Nibasis_event = length(aau(ibasis_event)),
    Nibasis_rcens = length(aau(ibasis_rcens))
  )

  #----- priors and hyperparameters

  # valid priors
  ok_dists <- nlist("normal",
                    student_t = "t",
                    "cauchy",
                    "hs",
                    "hs_plus",
                    "laplace",
                    "lasso") # disallow product normal
  ok_intercept_dists <- ok_dists[1:3]
  ok_aux_dists       <- ok_dists[1:3]
  ok_smooth_dists    <- c(ok_dists[1:3], "exponential")

  # priors
  user_prior_stuff <- prior_stuff <- lapply(seq_len(nt), function(k)
    handle_glm_prior(prior[[k]],
                   nvars = standata$s_K[k],
                   default_scale = 2.5,
                   link = NULL,
                   ok_dists = ok_dists)
  )
  user_prior_intercept_stuff <- prior_intercept_stuff <-
    lapply(seq_len(nt), function(k) handle_glm_prior(
      prior_intercept[[k]],
      nvars = 1,
      default_scale = 20,
      link = NULL,
      ok_dists = ok_intercept_dists)
    )
  user_prior_aux_stuff <- prior_aux_stuff <-
    lapply(seq_len(nt), function(k) handle_glm_prior(
      prior_aux[[k]],
      nvars = standata$s_vars[k],
      default_scale = get_default_aux_scale(basehaz[[k]]),
      link = NULL,
      ok_dists = ok_aux_dists)
    )


  # stop null priors if prior_PD is TRUE
  if (prior_PD) {
    if (is.null(prior))
      stop("'prior' cannot be NULL if 'prior_PD' is TRUE")
    if (is.null(prior_intercept && has_intercept) )
      stop("'prior_intercept' cannot be NULL if 'prior_PD' is TRUE")
    if (is.null(prior_aux))
      stop("'prior_aux' cannot be NULL if 'prior_PD' is TRUE")
  }

  # autoscaling of priors
  prior_stuff           <- lapply(seq_len(nt), function(k)
    autoscale_prior(prior_stuff[[k]], predictors = x[[k]]) )
  prior_intercept_stuff <- lapply(seq_len(nt), function(k)
    autoscale_prior(prior_intercept_stuff[[k]]))
  prior_aux_stuff       <- lapply(seq_len(nt), function(k)
    autoscale_prior(prior_aux_stuff[[k]]) )

  # priors
  standata$prior_dist              <- aa( sapply(prior_stuff, function(p) p$prior_dist) )
  standata$prior_dist_for_intercept <- aa( sapply(prior_intercept_stuff, function(i) i$prior_dist) )
  standata$prior_dist_for_aux      <-  aa( sapply(prior_aux_stuff, function(a) a$prior_dist ) )


  # hyperparameters
  standata$prior_mean    <- aau( lapply(prior_stuff, function(p) p$prior_mean ) )
  standata$prior_scale   <- aau( lapply(prior_stuff, function(p) p$prior_scale ) )
  standata$prior_df      <- aau( lapply(prior_stuff, function(p) p$prior_df ) )
  standata$prior_mean_for_intercept <- aau( lapply(prior_intercept_stuff, function(i) i$prior_mean) )
  standata$prior_scale_for_intercept <- aau( lapply(prior_intercept_stuff, function(i) i$prior_scale) )
  standata$prior_df_for_intercept   <- aau( lapply(prior_intercept_stuff, function(i) i$prior_df) )
  standata$prior_scale_for_aux <- aau( lapply( prior_aux_stuff, function(p) p$prior_scale))
  standata$prior_df_for_aux <- aau( lapply( prior_aux_stuff, function(p) p$prior_df))
  standata$global_prior_scale <- aau( lapply( prior_stuff, function(p) p$global_prior_scale))
  standata$global_prior_df <- aau( lapply( prior_stuff, function(p) p$global_prior_df))
  standata$slab_df <- aau( lapply( prior_stuff, function(p) p$slab_df))
  standata$slab_scale <- aau( lapply( prior_stuff, function(p) p$slab_scale))

  # standata$prior_mean_for_smooth01    <- prior_smooth_stuff01$prior_mean
  # standata$prior_scale_for_smooth01   <- prior_smooth_stuff01$prior_scale
  # standata$prior_df_for_smooth01      <- prior_smooth_stuff01$prior_df

  # any additional flags
  standata$prior_PD <- ai(prior_PD)

  #---------------
  # Prior summary
  #---------------

  prior_info <- lapply(seq_len(nt), function(k)
    summarize_jm_prior(
      user_priorEvent           = user_prior_stuff[[k]],
      user_priorEvent_intercept = user_prior_intercept_stuff[[k]],
      user_priorEvent_aux       = user_prior_aux_stuff[[k]],
      adjusted_priorEvent_scale           = prior_stuff[[k]]$prior_scale,
      adjusted_priorEvent_intercept_scale = prior_intercept_stuff[[k]]$prior_scale,
      adjusted_priorEvent_aux_scale       = prior_aux_stuff[[k]]$prior_scale,
      e_has_intercept  = standata$has_intercept[[k]],
      e_has_predictors = standata$s_K[[k]] > 0,
      basehaz = basehaz[[k]]
    )
    )

  #-----------
  # Fit model
  #-----------

  # obtain stan model code
  stanfit  <- stanmodels$mstte

  # specify parameters for stan to monitor
  stanpars <- c(if ( any(standata$has_intercept) ) "alpha",
                if ( any(standata$s_K) )            "beta",
                if ( any(standata$s_vars)  )       "aux")

  #
  # # fit model using stan
  if (algorithm == "sampling") { # mcmc
    args <- set_sampling_args(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      user_dots = list(...),
      user_adapt_delta = adapt_delta,
      show_messages = FALSE,
      prior = NULL)
    stanfit <- do.call(rstan::sampling, args)
  } else { # meanfield or fullrank vb
    args <- nlist(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      algorithm
    )
    args[names(dots)] <- dots
    stanfit <- do.call(rstan::vb, args)
  }

  check_stanfit(stanfit)


  # define new parameter names
  #  nms_tde    <- get_smooth_name(s_cpts, type = "smooth_coefs") # may be NULL
  #  nms_smooth <- get_smooth_name(s_cpts, type = "smooth_sd")    # may be NULL

  nms_int <- lapply(seq_len(nt), function(i) {
    if(has_intercept[[i]] > 0){
      append_trans(get_int_name_basehaz(basehaz[[i]]), i, transition_labels[i])
    } else {
      NULL
    }
  } )

  nms_beta <- lapply(seq_len(nt), function(i) {
    if(standata$s_K[i] > 0){
      append_trans(colnames(x[[i]]), i, transition_labels[i])
    } else {
      NULL
    }
  } )

  nms_aux <- lapply(seq_len(nt), function(i) {
    if(get_basehaz_name(basehaz[[i]]) != "exp"){
      append_trans(get_aux_name_basehaz(basehaz[[i]]), i, transition_labels[i])
    } else {
      NULL
    }
  } )

  nms_all <- ulist( c( nms_int, nms_beta, nms_aux) )

  nms_all <- c(nms_all, "log-posterior")
  # substitute new parameter names into 'stanfit' object
  stanfit <- replace_stanfit_nms(stanfit, nms_all)


  has_tde = rep(FALSE, nt)
  has_quadrature = rep( FALSE, nt) # not implemented

  # return an object of class 'stanmstte'
  fit <- nlist(stanfit,
               formula,
               has_tde = has_tde,
               has_quadrature = has_quadrature,
               data,
               transition_labels,
               model_frame      = mf,
               terms            = mt,
               xlevels          = lapply(seq_along(1:nt), function(i) .getXlevels(mt[[i]], mf[[i]]) ),
               x,
               s_cpts           = if (any(has_tde)) s_cpts else NULL,
               t_beg,
               t_end,
               status,
               event            = lapply(seq_along(1:nt), function(i) as.logical(status[[i]] == 1) ),
               delayed,
               basehaz,
               nobs             = lapply(seq_along(1:nt), function(i) nrow(mf[[i]]) ),
               nevents          = nevent,
               nlcens,
               nrcens,
               nicens,
               ncensor          = lapply(seq_along(1:nt), function(i) nlcens[[i]] + nrcens[[i]] + nicens[[i]]),
               ndelayed         = ndelay,
               prior_info,
               qnodes           = if (any(has_quadrature)) qnodes else NULL,
               algorithm,
               stan_function    = "mstte_stan",
               rstan_version    = utils::packageVersion("rstan"),
               call             = match.call(expand.dots = TRUE))

  msttestan(fit)
}


# Return the integer respresentation for the baseline hazard, used by Stan
#
# @param basehaz_name A character string, the type of baseline hazard.
# @return An integer, or NA if unmatched.
basehaz_for_stan <- function(basehaz_name) {
  switch(basehaz_name,
         weibull   = 1L,
         bs        = 2L,
         piecewise = 3L,
         ms        = 4L,
         exp       = 5L,
         gompertz  = 6L,
         NA)
}















