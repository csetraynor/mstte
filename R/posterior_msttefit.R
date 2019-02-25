#' Posterior predictions for muti-state time-to-event models
#' 
#' This function allows us to generate predicted quantities for survival
#' models at specified times. These quantities include the 
#' hazard rate, the cumulative hazard, or the survival probability.
#' Predictions are obtained using unique draws from the posterior distribution
#' of each of the model parameters and then summarised into a median and 
#' posterior uncertainty interval.
#' 
#' @export
#' @import splines2
#' 
#' 
#' @param newdataLong,newdataEvent Optionally, a data frame (or in the case of 
#'   \code{newdataLong} this can be a list of data frames) in which to look 
#'   for variables with which to predict. If omitted, the model matrices are used. 
#'   If new data is provided, then it should also contain the longitudinal 
#'   outcome data on which to condition when drawing the new group-specific 
#'   coefficients for individuals in the new data. Note that there is only
#'   allowed to be one row of data for each individual in \code{newdataEvent}, 
#'   that is, time-varying covariates are not allowed in the prediction data for
#'   the event submodel. Also, \code{newdataEvent} can optionally include a 
#'   variable with information about the last known survival time for the new
#'   individuals -- see the description for the \code{last_time} argument below
#'   -- however also note that when generating the survival probabilities it 
#'   is of course assumed that all individuals in \code{newdataEvent} have not 
#'   yet experienced the event (that is, any variable in \code{newdataEvent} that
#'   corresponds to the event indicator will be ignored).
#' @param type The type of prediction to return. The following are currently
#'   allowed:
#'   \itemize{
#'     \item \code{"surv"}: the estimated survival probability.
#'     \item \code{"cumhaz"}: the estimated cumulative hazard.
#'     \item \code{"haz"}: the estimated hazard rate.
#'   }
#' @param extrapolate A logical specifying whether to extrapolate the estimated 
#'   survival probabilities beyond the times specified in the \code{times} argument.
#'   If \code{TRUE} then the extrapolation can be further controlled using
#'   the \code{control} argument.
#' @param control A named list with parameters controlling extrapolation 
#'   of the estimated survival function when \code{extrapolate = TRUE}. The list
#'   can contain one or more of the following named elements: \cr
#'   \itemize{
#'     \item \code{epoints}: a positive integer specifying the number of  
#'     discrete time points at which to calculate the forecasted survival 
#'     probabilities. The default is 10.
#'     \item \code{edist}: a positive scalar specifying the amount of time 
#'     across which to forecast the estimated survival function, represented 
#'     in units of the time variable \code{time_var} (from fitting the model). 
#'     The default is to extrapolate between the times specified in the 
#'     \code{times} argument and the maximum event or censoring time in the 
#'     original data. If \code{edist} leads to times that are beyond
#'     the maximum event or censoring time in the original data then the 
#'     estimated survival probabilities will be truncated at that point, since
#'     the estimate for the baseline hazard is not available beyond that time.
#'   }
#' @param condition A logical specifying whether the estimated 
#'     subject-specific survival probabilities at time \code{t} should be 
#'     conditioned on survival up to a fixed time point \code{u}. The default 
#'     is for \code{condition} to be set to \code{TRUE}, unless standardised survival
#'     probabilities have been requested (by specifying \code{standardise = TRUE}),
#'     in which case \code{condition} must (and will) be set to \code{FALSE}.
#'     When conditional survival probabilities are requested, the fixed
#'     time point \code{u} will be either: (i) the value specified via the 
#'     \code{last_time} argument; or if the \code{last_time} argument is 
#'     \code{NULL} then the latest observation time for each individual 
#'     (taken to be the value in the \code{times} argument if \code{newdataEvent} 
#'     is specified, or the observed event or censoring time if \code{newdataEvent} 
#'     is \code{NULL}.
#' @param last_time A scalar, character string, or \code{NULL}. This argument 
#'     specifies the last known survival time for each individual when
#'     conditional predictions are being obtained. If 
#'     \code{newdataEvent} is provided and conditional survival predictions are being
#'     obtained, then the \code{last_time} argument can be one of the following:
#'     (i) a scalar, this will use the same last time for each individual in 
#'     \code{newdataEvent}; (ii) a character string, naming a column in 
#'     \code{newdataEvent} in which to look for the last time for each individual;
#'     (iii) \code{NULL}, in which case the default is to use the time of the latest 
#'     longitudinal observation in \code{newdataLong}. If \code{newdataEvent} is
#'     \code{NULL} then the \code{last_time} argument cannot be specified 
#'     directly; instead it will be set equal to the event or censoring time for
#'     each individual in the dataset that was used to estimate the model. 
#'     If standardised survival probabilities are requested (i.e. 
#'     \code{standardise = TRUE}) then conditional survival probabilities are
#'     not allowed and therefore the \code{last_time} argument is ignored.
#' @param ids An optional vector specifying a subset of IDs for whom the 
#'   predictions should be obtained. The default is to predict for all individuals
#'   who were used in estimating the model or, if \code{newdataLong} and 
#'   \code{newdataEvent} are specified, then all individuals contained in 
#'   the new data.
#' @param prob A scalar between 0 and 1 specifying the width to use for the 
#'   uncertainty interval (sometimes called credible interval) for the predictions. 
#'   For example \code{prob = 0.95} (the default) means that the 2.5th and 97.5th  
#'   percentiles will be provided.
#' @param times A scalar, a character string, or \code{NULL}. Specifies the  
#'   times at which the estimated survival probabilities should be calculated. 
#'   It can be either: (i) \code{NULL}, in which case it will default to the last known 
#'   survival time for each individual, as determined by the \code{last_time}
#'   argument; (ii) a scalar, specifying a time to estimate the survival probability
#'   for each of the individuals; or (iii) if \code{newdataEvent} is  
#'   provided, it can be the name of a variable in \code{newdataEvent} that 
#'   indicates the time at which the survival probabilities should be calculated  
#'   for each individual. 
#' @param standardise A logical specifying whether the estimated 
#'   subject-specific survival probabilities should be averaged
#'   across all individuals for whom the subject-specific predictions are 
#'   being obtained. This can be used to average over the covariate and random effects
#'   distributions of the individuals used in estimating the model, or the individuals 
#'   included in the \code{newdata} arguments. This approach of
#'   averaging across the observed distribution of the covariates is sometimes
#'   referred to as a "standardised" survival curve. If \code{standardise = TRUE}, 
#'   then the \code{times} argument must be specified and it must be constant across 
#'   individuals, that is, the survival probabilities must be calculated at the 
#'   same time for all individuals.
#' @param dynamic A logical that is only relevant if new data is provided
#'   via the \code{newdataLong} and \code{newdataEvent} arguments. If 
#'   \code{dynamic = TRUE}, then new group-specific parameters are drawn for 
#'   the individuals in the new data, conditional on their longitudinal 
#'   biomarker data contained in \code{newdataLong}. These group-specific
#'   parameters are then used to generate individual-specific survival probabilities
#'   for these individuals. These are often referred to as "dynamic predictions"
#'   in the joint modelling context, because the predictions can be updated
#'   each time additional longitudinal biomarker data is collected on the individual.
#'   On the other hand, if \code{dynamic = FALSE} then the survival probabilities
#'   will just be marginalised over the distribution of the group-specific
#'   coefficients; this will mean that the predictions will incorporate all
#'   uncertainty due to between-individual variation so there will likely be
#'   very wide credible intervals on the predicted survival probabilities.
#' @param scale A scalar, specifying how much to multiply the asymptotic 
#'   variance-covariance matrix for the random effects by, which is then
#'   used as the "width" (ie. variance-covariance matrix) of the multivariate
#'   Student-t proposal distribution in the Metropolis-Hastings algorithm. This
#'   is only relevant when \code{newdataEvent} is supplied and 
#'   \code{dynamic = TRUE}, in which case new random effects are simulated
#'   for the individuals in the new data using the Metropolis-Hastings algorithm.
#' @param draws An integer indicating the number of MCMC draws to return. 
#'   The default is to set the number of draws equal to 200, or equal to the 
#'   size of the posterior sample if that is less than 200. 
#' @param seed An optional \code{\link[=set.seed]{seed}} to use.
#' @param ... Currently unused.
#'
#' @details
#'   By default, the predicted quantities are evaluated conditional on observed 
#'   values of the fixed effect covariates. That is, predictions will be 
#'   obtained using either: 
#'   \itemize{
#'     \item the design matrices used in the original \code{\link{stan_surv}}
#'     or \code{\link{stan_jm}} model call, or
#'     \item the covariate values provided in the \code{newdata} argument
#'     (or \code{newdataLong} and \code{newdataEvent} arugments for the
#'     \code{stanjm} method). 
#'   }
#'   However, if you wish to average over the observed distribution 
#'   of the fixed effect covariates then this is possible -- such predictions
#'   are sometimes referred to as standardised survival probabilties -- see the 
#'   \code{standardise} argument.
#'   
#'   For \code{stansurv} objects, the predicted quantities are calculated for  
#'   each row of the prediction data, at the specified \code{times} as well as 
#'   any times generated through extrapolation (when \code{extrapolate = TRUE}).
#'   For \code{stanjm} objects, these quantities are calculated for each 
#'   individual, at the specified \code{times} as well as any times generated
#'   through extrapolation (when \code{extrapolate = TRUE}).
#'   
#'   The following also applies for \code{stanjm} objects.
#'   By default the survival probabilities are conditional on an individual's 
#'   group-specific coefficients (i.e. their individual-level random
#'   effects). If prediction data is provided via the \code{newdataLong}  
#'   and \code{newdataEvent} arguments, then the default behaviour is to
#'   sample new group-specific coefficients for the individuals in the  
#'   new data using a Monte Carlo scheme that conditions on their 
#'   longitudinal outcome data provided in \code{newdataLong} 
#'   (sometimes referred to as "dynamic predictions", see Rizopoulos
#'   (2011)). This default behaviour can be stopped by specifying 
#'   \code{dynamic = FALSE}, in which case the predicted survival
#'   probabilities will be marginalised over the distribution of the 
#'   group-specific coefficients. This has the benefit that the user does
#'   not need to provide longitudinal outcome measurements for the new 
#'   individuals, however, it does mean that the predictions will incorporate
#'   all the uncertainty associated with between-individual variation, since
#'   the predictions aren't conditional on any observed data for the individual.
#'   
#' @note 
#'   Note that if any variables were transformed (e.g. rescaled) in the data 
#'   used to fit the model, then these variables must also be transformed in 
#'   \code{newdataLong} and \code{newdataEvent}. This only applies if variables  
#'   were transformed before passing the data to one of the modeling functions and  
#'   \emph{not} if transformations were specified inside the model formula.
#'    
#' @return A data frame of class \code{survfit.stanjm}. The data frame includes 
#'   columns for each of the following: 
#'   (i) the median of the posterior predictions of the estimated survival
#'   probabilities (\code{survpred});
#'   (ii) each of the lower and upper limits of the corresponding uncertainty 
#'   interval for the estimated survival probabilities (\code{ci_lb} and 
#'   \code{ci_ub});
#'   (iii) a subject identifier (\code{id_var}), unless standardised survival
#'   probabilities were estimated;
#'   (iv) the time that the estimated survival probability is calculated for 
#'   (\code{time_var}).
#'   The returned object also includes a number of additional attributes.
#' 
#' @seealso 
#'   \code{\link{plot.survfit.stanjm}} for plotting the estimated survival  
#'   probabilities \cr
#'   \code{\link{ps_check}} for for graphical checks of the estimated 
#'   survival function \cr
#'   \code{\link{posterior_traj}} for estimating the
#'   marginal or subject-specific longitudinal trajectories \cr
#'   \code{\link{plot_stack_jm}} for combining plots of the estimated 
#'   subject-specific longitudinal trajectory and survival function
#'   
#' @references 
#'   Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in 
#'   joint models for longitudinal and time-to-event data. \emph{Biometrics}
#'   \strong{67}, 819.
#'      
#' @examples
#' \donttest{
#'   # Run example model if not already loaded
#'   if (!exists("example_jm")) example(example_jm)
#'   
#'   # Obtain subject-specific survival probabilities for a few
#'   # selected individuals in the estimation dataset who were  
#'   # known to survive up until their censoring time. By default
#'   # the posterior_survfit function will estimate the conditional
#'   # survival probabilities, that is, conditional on having survived
#'   # until the event or censoring time, and then by default will
#'   # extrapolate the survival predictions forward from there.  
#'   ps1 <- posterior_survfit(example_jm, ids = c(7,13,15))
#'   # We can plot the estimated survival probabilities using the
#'   # associated plot function
#'   plot(ps1)
#'   
#'   # If we wanted to estimate the survival probabilities for the
#'   # same three individuals as the previous example, but this time
#'   # we won't condition on them having survived up until their 
#'   # censoring time. Instead, we will estimate their probability
#'   # of having survived between 0 and 5 years given their covariates
#'   # and their estimated random effects.
#'   # The easiest way to achieve the time scale we want (ie, 0 to 5 years)
#'   # is to specify that we want the survival time estimated at time 0
#'   # and then extrapolated forward 5 years. We also specify that we
#'   # do not want to condition on their last known survival time.
#'   ps2 <- posterior_survfit(example_jm, ids = c(7,13,15), times = 0,
#'     extrapolate = TRUE, condition = FALSE, control = list(edist = 5))
#'     
#'   # Instead we may want to estimate subject-specific survival probabilities 
#'   # for a set of new individuals. To demonstrate this, we will simply take
#'   # the first two individuals in the estimation dataset, but pass their data
#'   # via the newdata arguments so that posterior_survfit will assume we are 
#'   # predicting survival for new individuals and draw new random effects 
#'   # under a Monte Carlo scheme (see Rizopoulos (2011)).
#'   ndL <- pbcLong[pbcLong$id %in% c(1,2),]
#'   ndE <- pbcSurv[pbcSurv$id %in% c(1,2),]
#'   ps3 <- posterior_survfit(example_jm,
#'     newdataLong = ndL, newdataEvent = ndE,
#'     last_time = "futimeYears", seed = 12345)
#'   head(ps3)
#'   # We can then compare the estimated random effects for these 
#'   # individuals based on the fitted model and the Monte Carlo scheme
#'   ranef(example_jm)$Long1$id[1:2,,drop=FALSE] # from fitted model
#'   colMeans(attr(ps3, "b_new"))                # from Monte Carlo scheme
#'   
#'   # Lastly, if we wanted to obtain "standardised" survival probabilities, 
#'   # (by averaging over the observed distribution of the fixed effect 
#'   # covariates, as well as averaging over the estimated random effects
#'   # for individuals in our estimation sample or new data) then we can 
#'   # specify 'standardise = TRUE'. We can then plot the resulting 
#'   # standardised survival curve.
#'   ps4 <- posterior_survfit(example_jm, standardise = TRUE, 
#'                            times = 0, extrapolate = TRUE)
#'   plot(ps4)
#' }
#'
posterior_msttefit <- function(object, ...) UseMethod("posterior_msttefit")

#' @rdname posterior_msttefit
#' @export
#'
posterior_msttefit.stanmstte <- function(object, 
                                       newdata     = NULL, 
                                       type        = "surv", 
                                       extrapolate = TRUE, 
                                       control     = list(), 
                                       condition   = NULL, 
                                       last_time   = NULL, 
                                       prob        = 0.95, 
                                       id_var      = NULL,
                                       times       = NULL,
                                       standardise = FALSE, 
                                       draws       = NULL, 
                                       seed        = NULL, 
                                       ...) {
  validate_stansurv_object(object)
  
  basehaz  <- object$basehaz
  
  if (!is.null(seed)) 
    set.seed(seed)
  if (requires.idvar(object) && is.null(id_var))
    STOP_id_var_required()
  
  dots <- list(...)
  
  newdata <- validate_newdata(newdata)
  has_newdata <- not.null(newdata)
  
  # Obtain a vector of unique subject ids 
  if (is.null(id_var)) {
    if (is.null(newdata)) {
      id_list <- seq(nrow(object$model_data))
    } else {
      id_list <- seq(nrow(newdata))
    }
  } else {
    if (is.null(newdata)) {
      id_list <- unique(object$model_data[[id_var]])
    } else {
      id_list <- unique(newdata[[id_var]])
    }
  }
  
  # Last known survival time for each individual
  if (is.null(newdata)) { # user did not specify newdata
    if (!is.null(last_time))
      stop("'last_time' cannot be provided when newdata is NULL, since times ",
           "are taken to be the event or censoring time for each individual.")
    last_time <- object$exittime
  } else { # user specified newdata
    if (is.null(last_time)) { # use latest longitudinal observation
      last_time <- rep(0, length(id_list))
    } else if (is.string(last_time)) {
      if (!last_time %in% colnames(ndE))
        stop("Cannot find 'last_time' column named in newdataEvent.")
      last_time <- newdata[[last_time]]      
    } else if (is.scalar(last_time)) {
      last_time <- rep(last_time, nrow(newdata)) 
    } else if (is.numeric(last_time) && (length(last_time) > 1L)) {
      last_time <- last_time[as.character(id_list)]
    } else {
      stop("Bug found: could not reconcile 'last_time' argument.")
    }
    names(last_time) <- as.character(id_list)
  }
  
  # Prediction times
  if (standardise) { # standardised survival probs
    times <- 
      if (is.null(times)) {
        stop("'times' cannot be NULL for obtaining standardised survival probabilities.")
      } else if (is.scalar(times)) {
        rep(times, length(id_list))
      } else {
        stop("'times' should be a numeric vector of length 1 in order to obtain ",
             "standardised survival probabilities (the subject-specific survival ",
             "probabilities will be calculated at the specified time point, and ",
             "then averaged).")      
      }    
  } else if (is.null(newdata)) { # subject-specific survival probs without newdata
    times <- 
      if (is.null(times)) {
        object$exittime
      } else if (is.scalar(times)) {
        rep(times, length(id_list))
      } else {
        stop("If newdata is NULL then 'times' must be NULL or a single number.")     
      }
  } else { # subject-specific survival probs with newdata
    times <- 
      if (is.null(times)) {
        times <- last_time
      } else if (is.scalar(times)) {
        rep(times, length(id_list))
      } else if (is.string(times)) {
        if (!times %in% colnames(ndE))
          stop("Variable specified in 'times' argument could not be found in newdata.")
        tapply(newdata[[times]], newdata[[id_var]], FUN = max)
      } else {
        stop("If newdata is specified then 'times' can only be the name of a ",
             "variable in newdata, or a single number.")
      }
  } 
  
  maxtime <- max(object$exittime)
  if (any(times > maxtime))
    stop("'times' are not allowed to be greater than the last event or censoring ",
         "time (since unable to extrapolate the baseline hazard).")
  
  
  # User specified extrapolation
  if (extrapolate) {
    control <- extrapolation_control(control, ok_args = c("epoints", "edist"))
    if (not.null(control$edist)) {
      endtime <- times + control$edist
    } else {
      endtime <- maxtime
    }
    endtime <- truncate(endtime, upper = maxtime)
    time_seq <- get_time_seq(control$epoints, times, endtime, simplify = FALSE)
  } else {
    time_seq <- list(times) # no extrapolation
  }
  
  # Conditional survival times
  if (is.null(condition)) {
    condition <- ifelse(type == "surv", !standardise, FALSE)
  } else if (condition && standardise) {
    stop("'condition' cannot be TRUE for standardised survival probabilities.")
  }
  
  # Get stanmat parameter matrix for specified number of draws
  stanmat <- sample_stanmat(object, draws = draws, default_draws = 200)
  pars    <- extract_pars(object, stanmat)
  
  # Calculate survival probability at each increment of extrapolation sequence
  surv <- lapply(time_seq, .pp_execute_surv, 
                 object      = object,
                 newdata     = newdata,
                 pars        = pars,
                 type        = type,
                 id_var      = id_var,
                 standardise = standardise)
  
  # If calculating log likelihood, then no further summarising required
  if (type == "ll")
    return(surv)
  
  # Calculate survival probability at last known survival time and then
  # use that to calculate conditional survival probabilities
  if (condition) {
    if (!type == "surv")
      stop("'condition' can only be set to TRUE for survival probabilities.")
    cond_surv <- .pp_execute_surv(last_time,
                                  object  = object,
                                  newdata = newdata,
                                  pars    = pars,
                                  type    = type,
                                  id_var  = id_var)
    surv <- lapply(surv, function(x) truncate(x / cond_surv, upper = 1))        
  } 
  
  # Summarise posterior draws to get median and CI
  if (is.null(id_var)) 
    id_var <- "id"
  out <- .pp_summarise_surv(surv        = surv,
                            prob        = prob,
                            type        = type,
                            id_var      = id_var,
                            standardise = standardise)
  
  # Add attributes
  structure(out,
            id_var      = attr(out, "id_var"),
            time_var    = attr(out, "time_var"),
            extrapolate = extrapolate, 
            control     = control, 
            condition   = condition, 
            standardise = standardise, 
            last_time   = last_time, 
            ids         = id_list, 
            draws       = draws, 
            seed        = seed, 
            class       = c("msttefit.stanmstte", "data.frame"))
}

#' @rdname posterior_msttefit
#' @export
#'
posterior_msttefit.stanmsjm <- function(object, newdataLong = NULL, newdataMs = NULL,
                                     extrapolate = TRUE, control = list(), 
                                     condition = NULL, last_time = NULL, prob = 0.95,  ids, times = NULL, standardise = FALSE, time_start_var = NULL,
                                     dynamic = TRUE, scale = 1.5,
                                     draws = NULL, seed = NULL, ...) {
  validate_stanmsjm_object(object)
  
  M        <- object$n_markers
  id_var   <- object$id_var
  time_var <- object$time_var
  basehaz  <- object$basehaz
  assoc    <- object$assoc
  family   <- object$family
  n_trans  <- object$n_trans
  transition_labels <- object$transition_labels
  obs_list <-  object$obs_list
  id_list <- object$id_list
  
  if (not.null(seed)) 
    set.seed(seed)
  if (missing(ids)) 
    ids <- NULL
  
  dots <- list(...)
  
  # Temporary stop, until make_assoc_terms can handle it
  sel_stop <- grep("^shared", rownames(object$assoc))
  if (any(unlist(object$assoc[sel_stop,])))
    stop("'posterior_survfit' cannot yet be used with shared_b or shared_coef ",
         "association structures.") 
  
  # Construct prediction data
  # ndL: dataLong to be used in predictions
  # ndE: dataEvent to be used in predictions
  if (!identical(is.null(newdataLong), is.null(newdataMs)))
    stop("Both newdataLong and newdataMs must be supplied together.")
  if (is.null(newdataLong)) { # user did not specify newdata
    dats <- get_model_data(object)
    ndL <- dats[1:M]
    ndE <- dats[["Event"]]
  } else { # user specified newdata
    if (!dynamic)
      stop2("Marginalised predictions for the event outcome are ",
            "not currently implemented.")
    has_newdata <- not.null(newdataMs)
    newdatas <- validate_newdatas_ms(object, newdataLong, newdataMs)
    ndL <- newdatas[1:M]
    ndMS <- newdatas[(M+1):(M+n_trans)]  
  }
  if(!is.null(newdataMs)){
    newdataMs_list <- lapply(seq_len(n_trans), function(i) newdataMs[newdataMs[["trans"]] == i, ])
    id_list <- lapply(newdataMs_list, function(x) extract_id(x, id_var))
    
    if(is.null(time_start_var) & !("Tstart" %in% colnames(newdataMs)) ){
      stop2("Provide a column name from newdataMS with time start variable.")
    } else if(is.null(time_start_var) ){
      warning2("Setting time_start_var to Tstart.")
      time_start_var <- "Tstart"
    }
    obs_list <-  mapply(FUN = match_obs,
                         newdataMs = newdataMs_list,
                         id_trans = id_list,
                         MoreArgs = nlist(long = newdataLong,
                                          time_var, id_var, time_start_var), SIMPLIFY = FALSE )
  }
  
  if (!is.null(ids)) { # user specified a subset of ids
    ndL <- subset_ids(object, ndL, ids)
    ndMS <- subset_ids(object, ndMS, ids)
  }  
  id_list <- lapply(id_list, function(d) factor(unique(d)) ) # order of ids from data, not ids arg
  
  # Last known survival time for each individual
  if (is.null(newdataLong)) { # user did not specify newdata
    if (!is.null(last_time))
      stop("'last_time' cannot be provided when newdata is NULL, since times ",
           "are taken to be the event or censoring time for each individual.")
    last_time <- object$eventtime[as.character(id_list)]
  } else { # user specified newdata
    if (is.null(last_time)) { # use latest longitudinal observation
      max_ytimes <- do.call("cbind", lapply(ndL, function(x) 
        tapply(x[[time_var]], x[[id_var]], FUN = max)))
      last_time <- apply(max_ytimes, 1L, max)
      # re-order last-time according to id_list
      last_time <- last_time[as.character(id_list)]
    } else if (is.character(last_time) && (length(last_time) == 1L)) {
      if (!last_time %in% colnames(ndMS))
        stop("Cannot find 'last_time' column named in newdataEvent.")
      last_time <- ndMS[[last_time]]      
    } else if (is.numeric(last_time) && (length(last_time) == 1L)) {
      last_time <- lapply(id_list, function(d) rep(last_time, length(d)) )
    } else if (is.numeric(last_time) && (length(last_time) > 1L)) {
      last_time <- last_time[as.character(id_list)]
    } else {
      stop("Bug found: could not reconcile 'last_time' argument.")
    }
    nms_lt <- lapply(id_list, function(d) as.character(d) )
    for(i in seq_len(n_trans)){
      names(last_time[[i]]) <- nms_lt[[i]]
    }
  }   
  
  # Prediction times
  if (standardise) { # standardised survival probs
    times <- 
      if (is.null(times)) {
        stop("'times' cannot be NULL for obtaining standardised survival probabilities.")
      } else if (is.numeric(times) && (length(times) == 1L)) {
        rep(times, length(id_list))
      } else {
        stop("'times' should be a numeric vector of length 1 in order to obtain ",
             "standardised survival probabilities (the subject-specific survival ",
             "probabilities will be calculated at the specified time point, and ",
             "then averaged).")      
      }    
  } else if (is.null(newdataLong)) { # subject-specific survival probs without newdata
    times <- 
      if (is.null(times)) {
        object$eventtime[as.character(id_list)]
      } else if (is.numeric(times) && (length(times) == 1L)) {
        rep(times, length(id_list))
      } else {
        stop("If newdata is NULL then 'times' must be NULL or a single number.")     
      }
  } else { # subject-specific survival probs with newdata
    times <- 
      if (is.null(times)) {
        times <- last_time
      } else if (is.character(times) && (length(times) == 1L)) {
        if (!times %in% colnames(ndMS))
          stop("Variable specified in 'times' argument could not be found in newdata.")
        tapply(ndMS[[times]], ndMS[[id_var]], FUN = max)
      } else if (is.numeric(times) && (length(times) == 1L)) {
        lapply(id_list, function(d) rep(times, length(d)) )
      } else {
        stop("If newdata is specified then 'times' can only be the name of a ",
             "variable in newdata, or a single number.")      
      }
  }
  if (!identical(length(times), length(id_list)))
    stop(paste0("length of the 'times' vector should be equal to the number of individuals ",
                "for whom predictions are being obtained (", length(id_list), ")."))     
  maxtime <- sapply(object$ms_mod, function(m) max(m$t_end) )
  lapply(seq_len(n_trans), function(k){
    if (any(times[[k]] > maxtime[k]))
      stop("'times' are not allowed to be greater than the last event or censoring ",
           "time (since unable to extrapolate the baseline hazard).")
  }
    )

  # User specified extrapolation
  if (extrapolate) {
    ok_args <- c("epoints", "edist")
    control <- extrapolation_control(control, ok_args = ok_args)
    endtime <- if (!is.null(control$edist)){
      lapply(times, function(t)  t + control$edist)
    } else {
      maxtime
      }
    for(k in n_trans){
      endtime[[k]][endtime[[k]] > maxtime[[k]]] <- maxtime[[k]] # nothing beyond end of baseline hazard 
    } 
    time_seq <- lapply(seq_len(n_trans), function(k){
      get_time_seq(control$epoints, times[[k]], endtime[[k]], simplify = FALSE)
    } )
  } else time_seq <- list(times) # no extrapolation
  
  # Conditional survival times
  if (is.null(condition)) {
    condition <- !standardise
  } else if (condition && standardise) {
    stop("'condition' cannot be set to TRUE if standardised survival ",
         "probabilities are requested.")
  }
  
  # Get stanmat parameter matrix for specified number of draws
  stanmat <- sample_stanmat(object, draws = draws, default_draws = 200)
  
  # Draw b pars for new individuals
  if (dynamic && !is.null(newdataMs)) {
    stanmat <- lapply(1:2, function(h) {
      simulate_b_pars_ms(object, stanmat = stanmat, ndL = ndL, ndE = ndMS[h],  ids = id_list[[h]], times = last_time[[h]], scale = scale) })
    b_new <- attr(stanmat, "b_new")
    acceptance_rate <- attr(stanmat, "acceptance_rate")
  }
  
  pars <- extract_pars(object, stanmat) # list of stanmat arrays
  
  # Matrix of surv probs at each increment of the extrapolation sequence
  # NB If no extrapolation then length(time_seq) == 1L
  pars_trans <- get_transitions_pars(pars, transition_labels)
  
  surv_t <- lapply(seq_len(n_trans), function(h) {
    lapply(time_seq[[h]], function(t) {  
    if (!identical(length(t), length(id_list[[h]])))
      stop("Bug found: the vector of prediction times is not the same length ",
           "as the number of individuals.")
    dat <- .pp_data_jm2(object, newdataLong = ndL[[1]][obs_list[[h]], ], newdataEvent = ndMS[[h]],  ids = id_list[[h]], etimes = t, long_parts = FALSE)
    surv_t <- .ll_survival(object, data = dat, pars = pars[[h]], survprob = TRUE)
    if (is.vector(surv_t) == 1L) 
      surv_t <- t(surv_t)   # transform if only one individual
    surv_t[, (t == 0)] <- 1 # avoids possible NaN due to numerical inaccuracies
    if (standardise) {      # standardised survival probs
      surv_t <- matrix(rowMeans(surv_t), ncol = 1)
      dimnames(surv_t) <- list(iterations = NULL, "standardised_survprob") 
    } else {
      dimnames(surv_t) <- list(iterations = NULL, ids = id_list[[h]])
    }
    surv_t
  })
  })
  
  # If conditioning, need to obtain matrix of surv probs at last known surv time
  if (condition) {
    cond_dat <- .pp_data_jm(object, newdataLong = ndL, newdataEvent = ndMS, 
                            ids = id_list, etimes = last_time, long_parts = FALSE)
    # matrix of survival probs at last_time 
    cond_surv <- .ll_survival(object, data = cond_dat, pars = pars, survprob = TRUE)
    if (is.vector(cond_surv) == 1L)
      cond_surv <- t(cond_surv)        # transform if only one individual
    cond_surv[, (last_time == 0)] <- 1 # avoids possible NaN due to numerical inaccuracies
    surv <- lapply(surv_t, function(x) { # conditional survival probs
      vec <- x / cond_surv
      vec[vec > 1] <- 1 # if t was before last_time then surv prob may be > 1
      vec
    })        
  } else surv <- surv_t
  
  # Summarise posterior draws to get median and ci
  out <- do.call("rbind", lapply(
    seq_along(surv), function(x, standardise, id_list, time_seq, prob) {
      val <- median_and_bounds(surv[[x]], prob, na.rm = TRUE)
      if (standardise) {
        data.frame(TIMEVAR = unique(time_seq[[x]]), val$med, val$lb, val$ub)        
      } else
        data.frame(IDVAR = id_list, TIMEVAR = time_seq[[x]], val$med, val$lb, val$ub) 
    }, standardise, id_list, time_seq, prob))
  out <- data.frame(out)
  colnames(out) <- c(if ("IDVAR" %in% colnames(out)) id_var,
                     time_var, "survpred", "ci_lb", "ci_ub")  
  if (id_var %in% colnames(out)) { # data has id column -- sort by id and time
    out <- out[order(out[, id_var, drop = F], out[, time_var, drop = F]), , drop = F]
  } else { # data does not have id column -- sort by time only
    out <- out[order(out[, time_var, drop = F]), , drop = F]
  }
  rownames(out) <- NULL
  
  # temporary hack so that predictive_error can call posterior_survfit
  # with two separate conditioning times...
  fn <- tryCatch(sys.call(-1)[[1]], error = function(e) NULL)
  if (!is.null(fn) && 
      grepl("predictive_error", deparse(fn), fixed = TRUE) &&
      "last_time2" %in% names(dots)) {
    last_time2 <- ndMS[[dots$last_time2]]
    cond_dat2 <- .pp_data_jm(object, newdataLong = ndL, newdataEvent = ndMS, 
                             ids = id_list, etimes = last_time2, long_parts = FALSE)
    cond_surv2 <- .ll_survival(object, data = cond_dat2, pars = pars, survprob = TRUE)
    if (is.vector(cond_surv2) == 1L)
      cond_surv2 <- t(cond_surv2)        # transform if only one individual
    cond_surv2[, (last_time2 == 0)] <- 1 # avoids possible NaN due to numerical inaccuracies
    surv2 <- lapply(surv_t, function(x) { # conditional survival probs
      vec <- x / cond_surv2
      vec[vec > 1] <- 1 # if t was before last_time then surv prob may be > 1
      vec
    })
    out2 <- do.call("rbind", lapply(
      seq_along(surv2), function(x, standardise, id_list, time_seq, prob) {
        val <- median_and_bounds(surv2[[x]], prob, na.rm = TRUE)
        data.frame(IDVAR = id_list, TIMEVAR = time_seq[[x]], val$med) 
      }, standardise, id_list, time_seq, prob))
    out2 <- data.frame(out2)
    colnames(out2) <- c(id_var, time_var, "survpred_eventtime")  
    out2 <- out2[order(out2[, id_var, drop = F], out2[, time_var, drop = F]), , drop = F]
    rownames(out2) <- NULL
    out <- merge(out, out2)
  }
  
  class(out) <- c("msttefit.stanmsjm", "data.frame")
  out <- structure(out, id_var = id_var, time_var = time_var, extrapolate = extrapolate, 
                   control = control, standardise = standardise, condition = condition, 
                   last_time = last_time, ids = id_list, draws = draws, seed = seed, 
                   offset = offset)
  if (dynamic && has_newdata) {
    out <- structure(out, b_new = b_new, acceptance_rate = acceptance_rate)
  }
  out
}



# internal ----------------------------------------------------------------

# Return a list with the control arguments for interpolation and/or
# extrapolation in posterior_predict.stanmvreg and posterior_survfit.stanjm
#
# @param control A named list, being the user input to the control argument
#   in the posterior_predict.stanmvreg or posterior_survfit.stanjm call
# @param ok_args A character vector of allowed control arguments
# @return A named list
extrapolation_control <- 
  function(control = list(), ok_args = c("epoints", "edist", "eprop")) {
    defaults <- list(ipoints = 15, epoints = 15, edist = NULL, eprop = 0.2, last_time = NULL)
    if (!is.list(control)) {
      stop("'control' should be a named list.")
    } else if (!length(control)) {
      control <- defaults[ok_args] 
    } else {  # user specified control list
      nms <- names(control)
      if (!length(nms))
        stop("'control' should be a named list.")
      if (any(!nms %in% ok_args))
        stop(paste0("'control' list can only contain the following named arguments: ",
                    paste(ok_args, collapse = ", ")))
      if (all(c("edist", "eprop") %in% nms))
        stop("'control' list cannot include both 'edist' and 'eprop'.")        
      if (("ipoints" %in% ok_args) && is.null(control$ipoints)) 
        control$ipoints <- defaults$ipoints   
      if (("epoints" %in% ok_args) && is.null(control$epoints)) 
        control$epoints <- defaults$epoints  
      if (is.null(control$edist) && is.null(control$eprop)) 
        control$eprop <- defaults$eprop
    }
    return(control)
  }
