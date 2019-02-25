# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017 Sam Brilleman
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Posterior predictions for survival models
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
#' @templateVar stanregArg object#'
#' @param newdata Optionally, a data frame in which to look for variables with
#'   which to predict. If omitted, the model matrix is used. If \code{newdata}
#'   is provided and any variables were transformed (e.g. rescaled) in the data
#'   used to fit the model, then these variables must also be transformed in
#'   \code{newdata}. This only applies if variables were transformed before
#'   passing the data to one of the modeling functions and \emph{not} if
#'   transformations were specified inside the model formula. Also,
#'   \code{newdata} can optionally include a variable with information
#'   about the last known survival time for the new individuals --
#'   see the description for the \code{last_time} argument below
#'   -- however also note that when generating the survival probabilities it
#'   is of course assumed that all individuals in \code{newdata} have not
#'   yet experienced the event (that is, any variable in \code{newdataEvent}
#'   that corresponds to the event indicator will be ignored).
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
#'   were transformed before passing the data to one of the modeling functions
#'   and \emph{not} if transformations were specified inside the model formula.
#'
#' @return A data frame of class \code{survfit.stanjm}. The data frame includes
#'   columns for each of the following:
#'   (i) the median of the posterior predictions (\code{median});
#'   (ii) each of the lower and upper limits of the corresponding uncertainty
#'   interval for the posterior predictions (\code{ci_lb} and \code{ci_ub});
#'   (iii) an observation identifier (for \code{stan_surv} models) or an
#'   individual identifier (for \code{stan_jm} models), unless standardised
#'   predictions were requested;
#'   (iv) the time that the prediction corresponds to.
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
posterior_mstte_fit <- function(object, ...) UseMethod("posterior_mstte_fit")


#' @rdname posterior_mstte_fit
#' @method posterior_mstte_fit stanmstte
#' @export
#'
posterior_mstte_fit.stanmstte <- function(object,
                                  newdata     = NULL,
                                  type        = "surv",
                                  extrapolate = TRUE,
                                  control     = list(),
                                  condition   = NULL,
                                  last_time   = NULL,
                                  prob        = 0.95,
                                  times       = NULL,
                                  standardise = FALSE,
                                  draws       = NULL,
                                  seed        = NULL,
                                  ids = NULL,
                                  ...) {

  validate_stanmstte_object(object)
  if (!is.null(seed))
    set.seed(seed)

  N = max(newdata$trans)

  labs = object$transition_labels
  last_time <- maybe_broadcast(last_time, N)
  times <- maybe_broadcast(times, N)
  condition <- maybe_broadcast(condition, N)

  newdata = lapply(seq_len(N), function(i)
    newdata[newdata$trans == i, ] )

  dots <- list(...)

  # Get stanmat parameter matrix for specified number of draws
  stanmat <- sample_stanmat(object, draws = draws, default_draws = NA)

  out <- list()

  for(h in seq_len(N)){

    basehaz  <- object$basehaz[[h]]
    newdata_h <- newdata[[h]]
    times_h <- times[[h]]
    last_time_h <- last_time[[h]]

    if (is.null(newdata_h) && object$ndelayed[[h]])
      stop("Prediction data for 'posterior_survfit' cannot include delayed ",
           "entry. If you estimated a model with delayed entry, you will ",
           "not be able to obtain predictions using the estimation data frame. ",
           "You must provide prediction data via the 'newdata' argument, and ",
           "indicate delayed entry via the 'last_time' argument.")

    newdata_h <- validate_newdata(newdata_h)
    has_newdata <- not.null(newdata_h)

    # Obtain a vector of unique subject ids
    if (is.null(newdata_h)) {
      id_list <- seq(nrow(get_model_data(object)[[h]] ))
    } else {
      id_list <- seq(nrow(newdata_h))
    }

    # Last known survival time for each individual
    if (is.null(newdata_h)) { # user did not specify newdata
      if (!is.null(last_time_h))
        stop("'last_time' cannot be provided when newdata is NULL, since times ",
             "are taken to be the event or censoring time for each individual.")
      last_time_h <- object$eventtime[[h]]
    } else { # user specified newdata
      if (is.null(last_time_h)) { # assume at risk from time zero
        last_time_h <- rep(0, length(id_list))
      } else if (is.string(last_time_h)) {
        if (!last_time_h %in% colnames(newdata_h))
          stop("Cannot find 'last_time' column named in newdata")
        last_time_h <- newdata_h[[last_time_h]]
      } else if (is.scalar(last_time_h)) {
        last_time_h <- rep(last_time_h, nrow(newdata_h))
      } else if (any(!is.numeric(last_time_h), !length(last_time_h) == nrow(newdata_h))) {
        stop("Bug found: could not reconcile 'last_time' argument.")
      }
      names(last_time_h) <- as.character(id_list)
    }

    # Prediction times
    if (standardise) { # standardised survival probs
      times_h <-
        if (is.null(times_h)) {
          stop("'times' cannot be NULL for obtaining standardised survival probabilities.")
        } else if (is.scalar(times_h)) {
          rep(times_h, length(id_list))
        } else {
          stop("'times' should be a numeric vector of length 1 in order to obtain ",
               "standardised survival probabilities (the subject-specific survival ",
               "probabilities will be calculated at the specified time point, and ",
               "then averaged).")
        }
    } else if (is.null(newdata_h)) { # subject-specific survival probs without newdata
      times_h <-
        if (is.null(times_h)) {
          object$eventtime[[h]]
        } else if (is.scalar(times_h)) {
          rep(times_h, length(id_list))
        } else {
          stop("If newdata is NULL then 'times' must be NULL or a single number.")
        }
    } else { # subject-specific survival probs with newdata
      times <-
        if (is.null(times_h)) {
          times_h <- last_time_h
        } else if (is.scalar(times_h)) {
          rep(times_h, length(id_list))
        } else if (is.string(times_h)) {
          if (!times_h %in% colnames(newdata_h))
            stop("Variable specified in 'times' argument could not be found in newdata.")
          times_h <- newdata_h[[times_h]]
        } else {
          stop("If newdata is specified then 'times' can only be the name of a ",
               "variable in newdata, or a single number.")
        }
    }
    maxtime <- max(object$eventtime[[h]])
    if (any(times_h > maxtime))
      stop("'times' are not allowed to be greater than the last event or ",
           "censoring time (since unable to extrapolate the baseline hazard).")

    # User specified extrapolation
    if (extrapolate) {
      control_h <- extrapolation_control(control, ok_args = c("epoints", "edist"))
      if (not.null(control_h$edist)) {
        endtime <- times_h + control_h$edist
      } else {
        endtime <- maxtime
      }
      endtime <- truncate(endtime, upper = maxtime)
      time_seq <- get_time_seq(control_h$epoints, times_h, endtime, simplify = FALSE)
    } else {
      time_seq <- list(times_h) # no extrapolation
    }


    # Conditional survival times
    if (is.null(condition[[h]])) {
      condition[[h]] <- ifelse(type == "surv", !standardise, FALSE)
    } else if (condition[[h]] && standardise) {
      stop("'condition' cannot be TRUE for standardised survival probabilities.")
    }

    # Get stanmat parameter matrix for transition
    string_to_grep <- paste0("*.trans\\(", labs[h], "\\)$")
    stanmat_h <- stanmat[ , grepl(string_to_grep, colnames(stanmat)), drop = FALSE]
    pars    <- extract_pars(object, stanmat)
    alpha <- pars$alpha[ , grepl(string_to_grep, colnames(pars$alpha)), drop = FALSE]
    beta <- pars$beta[ , grepl(string_to_grep, colnames(pars$beta)), drop = FALSE]
    beta_tde <- pars$beta_tde[ , grepl(string_to_grep, colnames(pars$beta_tde)), drop = FALSE]
    aux <- pars$aux[ , grepl(string_to_grep, colnames(pars$aux)), drop = FALSE]
    smooth <-pars$smooth[ , grepl(string_to_grep, colnames(pars$smoot)), drop = FALSE]

    pars_h <- nlist(alpha, beta, beta_tde, aux, smooth, stanmat = stanmat_h)


    # Calculate survival probability at each increment of extrapolation sequence
    surv <- lapply(time_seq, .pp_calculate_surv,
                   object      = object,
                   newdata     = newdata_h,
                   pars        = pars_h,
                   type        = type,
                   standardise = standardise,
                   h = h)


    # Calculate survival probability at last known survival time and then
    # use that to calculate conditional survival probabilities
    if (condition[[h]]) {
      if (!type == "surv")
        stop("'condition' can only be set to TRUE for survival probabilities.")
      cond_surv <- .pp_calculate_surv(last_time_h,
                                      object  = object,
                                      newdata = newdata_h,
                                      pars    = pars_h,
                                      type    = type,
                                      h = h)
      surv <- lapply(surv, function(x) truncate(x / cond_surv, upper = 1))
    }

    # Summarise posterior draws to get median and CI
    out[[h]] <- .pp_summarise_surv(surv        = surv,
                                   prob        = prob,
                                   standardise = standardise)

    # Add attributes
    out[[h]] <- structure(out[[h]],
                          id_var      = attr(out[[h]], "id_var"),
                          time_var    = attr(out[[h]], "time_var"),
                          type        = type,
                          extrapolate = extrapolate,
                          control     = control,
                          condition   = condition,
                          standardise = standardise,
                          last_time   = last_time,
                          ids         = id_list,
                          draws       = draws,
                          seed        = seed,
                          class       = c("survfit.stansurv", "data.frame"))

  }
  structure(out,
            seed        = seed,
            class       = c("posterior_mstte_fit.stanmstte", "list" )
  )
}

# -----------------  internal  ------------------------------------------------

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


# Calculate the desired prediction (e.g. hazard, cumulative hazard, survival
# probability) at the specified times
.pp_calculate_surv <- function(times,
                               object,
                               newdata      = NULL,
                               newdataLong  = NULL,
                               newdataEvent = NULL,
                               pars,
                               type         = "surv",
                               id_list      = NULL,
                               standardise  = FALSE,
                               h = h) {

  # Determine whether prediction type requires quadrature
  needs_quadrature <- type %in% c("cumhaz",
                                  "surv",
                                  "cdf",
                                  "logcumhaz",
                                  "logsurv",
                                  "logcdf")

  # Evaluate hazard, cumulative hazard, survival or failure probability
  ppdat <- .pp_data_mstte(object,
                        newdata       = newdata,
                        times         = times,
                        at_quadpoints = needs_quadrature,
                        k = h)

  out <- .pp_predict_surv(object,
                          data = ppdat,
                          pars = pars,
                          type = type,
                          h = h)

  # Transform if only one individual
  out <- transpose_vector(out)

  # Set survival probability == 1 if time == 0 (avoids possible NaN)
  if (type == "surv")
    out <- replace_where(out, times == 0, replacement = 1, margin = 2L)

  # Standardisation: within each iteration, calculate mean across individuals
  if (standardise) {
    out   <- row_means(out)
    ids   <- "standardised_survprob"
    times <- unique(times)
  } else {
    ids   <- if (is.null(id_list)) seq(ncol(out)) else id_list
  }
  dimnames(out) <- list(iterations = NULL, ids = ids)

  # Add subject ids and prediction times as an attribute
  structure(out, ids = ids, times = times)
}


# Evaluate hazard, cumulative hazard, survival or failure probability
#
# @param object A stansurv or stanjm object.
# @param data Output from .pp_data_surv or .pp_data_jm.
# @param pars Output from extract_pars.
# @param type The type of prediction quantity to return.
.pp_predict_surv <- function(object, ...) UseMethod(".pp_predict_surv")


.pp_predict_surv.stanmstte <- function(object,
                                     data,
                                     pars,
                                     type = "surv",
                                     h) {

  args <- nlist(basehaz   = object$basehaz[[h]],
                intercept = pars$alpha,
                betas     = pars$beta,
                betas_tde = pars$beta_tde,
                aux       = pars$aux,
                times     = data$pts,
                x         = data$x,
                s         = data$s)

  if (type %in% c("loghaz", "haz")) {
    # evaluate hazard; quadrature not relevant
    lhaz <- do.call(evaluate_log_haz, args)
  } else if (!data$has_quadrature[h]){
    # evaluate survival; without quadrature
    lsurv <- do.call(evaluate_log_surv, args)
  } else {
    # evaluate survival; with quadrature
    lhaz  <- do.call(evaluate_log_haz, args)
    lsurv <- -quadrature_sum(exp(lhaz), qnodes = data$qnodes, qwts = data$wts)
  }

  switch(type,
         loghaz    = lhaz,
         logcumhaz = log(-lsurv),
         logsurv   = lsurv,
         logcdf    = log(1 - exp(lsurv)),
         haz       = exp(lhaz),
         surv      = exp(lsurv),
         cumhaz    = -lsurv,
         cdf       = 1 - exp(lsurv),
         stop("Invalid input to the 'type' argument."))
}



# Summarise predictions into median, lower CI, upper CI
#
# @details Convert a list of matrices (with each element being a S by N matrix,
#   where S is the number of MCMC draws and N the number of individuals)
#   and collapse it across the MCMC iterations summarising it into median
#   and CI. The result is a data frame with K times N rows, where K was
#   the length of the original list.
.pp_summarise_surv <- function(surv,
                               prob        = NULL,
                               id_var      = NULL,
                               time_var    = NULL,
                               standardise = FALSE,
                               colnames    = NULL) {

  # Default variable names if not provided by the user
  if (is.null(id_var))
    id_var <- "id"
  if (is.null(time_var))
    time_var <- "time"

  # Extract ids and times for the predictions
  ids   <- uapply(surv, attr, "ids")
  times <- uapply(surv, attr, "times")

  # Determine the quantiles corresponding to the median and CI limits
  if (is.null(prob)) {
    probs <- 0.5 # median only
    nms   <- c(id_var, time_var, "median")
  } else {
    probs <- c(0.5, (1 - prob)/2, (1 + prob)/2) # median and CI
    nms   <- c(id_var, time_var, "median", "ci_lb", "ci_ub")
  }

  # Possibly overide default variable names for the returned data frame
  if (!is.null(colnames)) {
    nms <- c(id_var, time_var, colnames)
  }

  # Calculate mean and CI at each prediction time
  out <- data.frame(do.call("rbind", lapply(surv, col_quantiles_, probs)))
  out <- mutate_(out, id_var = ids, time_var = times)
  out <- row_sort(out, id_var, time_var)
  out <- col_sort(out, id_var, time_var)
  out <- set_rownames(out, NULL)
  out <- set_colnames(out, nms)

  # Drop excess info if standardised predictions were calculated
  if (standardise) { out[[id_var]] <- NULL; id_var <- NULL }

  structure(out,
            id_var   = id_var,
            time_var = time_var)
}


# ------------  print methods  ------------------------------------------------

print.survfit.stansurv <- function(x, digits = 4, ...) {

  x <- as.data.frame(x)
  sel <- c(attr(x, "time_var"), "median", "ci_lb", "ci_ub")
  for (i in sel)
    x[[i]] <- format(round(x[[i]], digits), nsmall = digits)


  cat(" num. individuals:", length(attr(x, "ids")), "\n")
  cat(" prediction type: ", tolower(get_survpred_name(attr(x, "type"))), "\n")
  cat(" standardised?:   ", yes_no_string(attr(x, "standardise")), "\n\n")
  print(x, quote = FALSE)
  invisible(x)
}


print.survfit.stanjm <- function(x, digits = 4, ...) {

  x <- as.data.frame(x)
  sel <- c(attr(x, "time_var"), "median", "ci_lb", "ci_ub")
  for (i in sel)
    x[[i]] <- format(round(x[[i]], digits), nsmall = digits)

  cat("stan_jm predictions\n")
  cat(" num. individuals:", length(attr(x, "ids")), "\n")
  cat(" prediction type: ", tolower(get_survpred_name(attr(x, "type"))), "\n")
  cat(" standardised?:   ", yes_no_string(attr(x, "standardise")), "\n\n")
  print(x, quote = FALSE)
  invisible(x)
}


print.fit.stanidm <- function(x, labels) {
  if(missing(labels)){
    labels <- handle_labels(x)
  }

  cat("stan_idm predictions\n")
  for(h in seq_along(x)){
    cat("\n ----------------- \n")
    cat(" ", labels[h], "\n" )
    print.survfit.stansurv(x[[h]])
  }
}

#' Method for stanmstte
#' @method as.data.frame posterior_mstte_fit.stanmstte
#' @export
as.data.frame.posterior_mstte_fit.stanmstte <- function(x, ...) {
  x <-
    mapply(`[<-`, x, 'transition', value = seq_along(x), SIMPLIFY = FALSE)

  do.call(rbind.data.frame, unclass(x), ...)
}


# -----------------  plot methods  --------------------------------------------

plot.fit.stanidm <- function(x,
                             ids    = NULL,
                             limits = c("ci", "none"),
                             xlab   = NULL,
                             ylab   = NULL,
                             facet_scales = "free",
                             ci_geom_args = NULL,
                             labels = "auto",
                             ...) {

  limits <- match.arg (limits)
  ci     <- as.logical(limits == "ci")

  if(is.null(xlab)) xlab <- lapply(seq_along(x), function(x) NULL)
  if(is.null(ylab)) ylab <- lapply(seq_along(x), function(x) NULL)
  if(is.null(ids)) ids <- lapply(seq_along(x), function(x) NULL)


  plotlist <- lapply(seq_along(x), function(i){

    x_n <- x[[i]]

    type        <- attr(x[[i]], "type")
    standardise <- attr(x[[i]], "standardise")
    id_var      <- attr(x[[i]], "id_var")
    time_var    <- attr(x[[i]], "time_var")

    if (is.null(xlab[[i]])) xlab[[i]] <- paste0("Time (", time_var, ")")
    if (is.null(ylab[[i]])) ylab[[i]] <- get_survpred_name(type)

    if (!is.null(ids[[i]])) {
      if (standardise)
        stop("'ids' argument cannot be specified when plotting standardised ",
             "survival probabilities.")
      x_n <- subset_ids(x[[i]], ids[[i]], id_var)
    } else {
      ids <- if (!standardise) attr(x_n, "ids") else NULL
    }
    if (!standardise) x_n$id <- factor(x_n[[id_var]])
    x_n$time <- x_n[[time_var]]

    geom_defaults <- list(color = "black")
    geom_mapp     <- list(mapping = aes_string(x = "time",
                                               y = "median"))
    geom_args     <- do.call("set_geom_args",
                             c(defaults = list(geom_defaults), list(...)))

    lim_defaults  <- list(alpha = 0.3)
    lim_mapp      <- list(mapping = aes_string(x = "time",
                                               ymin = "ci_lb",
                                               ymax = "ci_ub"))
    lim_args      <- do.call("set_geom_args",
                             c(defaults = list(lim_defaults), ci_geom_args))

    if ((!standardise) && (length(ids) > 60L))
      stop2("Too many individuals to plot for. Perhaps consider limiting ",
            "the number of individuals by specifying the 'ids' argument.")

    graph_base <-
      ggplot(x_n) +
      theme_bw() +
      coord_cartesian(ylim = get_survpred_ylim(type)) +
      do.call("geom_line", c(geom_mapp, geom_args))

    graph_facet <-
      if ((!standardise) && (length(ids) > 1L)) {
        facet_wrap(~ id, scales = facet_scales)
      } else NULL

    graph_limits <-
      if (ci) {
        do.call("geom_ribbon", c(lim_mapp, lim_args))
      } else NULL

    graph_labels <- labs(x = xlab[[i]], y = ylab[[i]])

    gg        <- graph_base + graph_facet + graph_limits + graph_labels
    class_gg  <- class(gg)
    class(gg) <- c("plot.survfit.stanjm", class_gg)
    gg
  })
  cowplot::plot_grid(
    plotlist = plotlist,
    axis = "rlbt",
    labels = labels,
    label_size = 10)

}


plot_stack_jm <- function(yplot, survplot) {

  if (!is(yplot, "list")) yplot <- list(yplot)

  lapply(yplot, function(x) {
    if (!is(x, "plot.predict.stanjm"))
      stop("'yplot' should be an object of class 'plot.predict.stanjm', ",
           "or a list of such objects.", call. = FALSE)
  })
  if (!is(survplot, "plot.survfit.stanjm"))
    stop("'survplot' should be an object of class 'plot.survfit.stanjm'.",
         call. = FALSE)

  y_build <- lapply(yplot, ggplot_build)
  y_layout <- lapply(y_build, function(x) x$layout$panel_layout)
  y_ids <- lapply(y_layout, function(x)
    if (!"id" %in% colnames(x)) NULL else x[["id"]])

  e_build <- ggplot_build(survplot)
  e_layout <- e_build$layout$panel_layout
  e_ids <- if (!"id" %in% colnames(e_layout)) NULL else e_layout[["id"]]

  if (!is.null(e_ids)) {
    lapply(y_ids, function(x, e_ids) {
      if (!all(sort(x) == sort(e_ids))) {
        stop("The individuals in the 'yplot' and 'survplot' appear to differ. Please ",
             "reestimate the plots using a common 'ids' argument.", call. = FALSE)
      }
    }, e_ids = e_ids)
  }

  vline <- lapply(seq_along(y_build), function(m) {
    L <- length(y_build[[m]]$data)
    dat <- y_build[[m]]$data[[L]]
    if (!"xintercept" %in% colnames(dat)) {
      found <- FALSE
    } else {
      found <- TRUE
      dat <- dat[, c("PANEL", "xintercept"), drop = FALSE]
      if (NROW(y_layout[[m]]) > 1) {
        panel_id_map <- y_layout[[m]][, c("PANEL", "id"), drop = FALSE]
        dat <- merge(dat, panel_id_map, by = "PANEL")
      }
      dat <- dat[, grep("PANEL", colnames(dat), invert = TRUE), drop = FALSE]
      colnames(dat) <- gsub("xintercept", paste0("xintercept", m), colnames(dat), fixed = TRUE)
    }
    list(dat = dat, found = found)
  })
  vline_found <- any(sapply(vline, function(x) x$found))
  if (!vline_found)
    cat("Could not find vertical line indicating last observation time in the",
        "plot of the longitudinal trajectory; you may wish to plot the longitudinal",
        "trajectories again with 'vline = TRUE' to aid interpretation.")
  vline_dat <- lapply(vline, function(x) x$dat)
  vline_alldat <- Reduce(function(...) merge(..., all = TRUE), vline_dat)
  vline_alldat$xintercept_max <-
    apply(vline_alldat[, grep("id", colnames(vline_alldat), invert = TRUE), drop = FALSE], 1, max)

  xmax <- max(sapply(c(y_build, list(e_build)), function(i) max(i$data[[1]]$x)))

  if ((!is.null(e_ids)) && (length(e_ids) > 20L)) {
    stop("Unable to generate 'plot_stack_jm' for this many individuals.", call. = FALSE)
  } else if ((!is.null(e_ids)) && (length(e_ids) > 3L)) {
    warning("'plot_stack_jm' is unlikely to be legible with more than a few individuals.",
            immediate. = TRUE, call. = FALSE)
  }

  if (!is.null(e_ids)) {
    graph_facet <- facet_wrap(~ id, scales = "free", nrow = 1)
  } else {
    graph_facet <- NULL
  }

  if (vline_found) {
    graph_vline <- geom_vline(aes_string(xintercept = "xintercept_max"),
                              vline_alldat, linetype = 2)
  } else {
    graph_vline <- NULL
  }

  graph_xlims <- expand_limits(x = c(0, xmax))

  survplot_updated <- survplot + graph_xlims + graph_facet + graph_vline

  yplot_updated <- lapply(yplot, function(x) x + graph_xlims + graph_facet)

  bayesplot::bayesplot_grid(
    plots = c(yplot_updated, list(survplot_updated)),
    grid_args = list(ncol = 1)
  )
}


# -----------------  helpers  ------------------------------------------------

# Return a user-friendly name for the prediction type
get_survpred_name <- function(x) {
  switch(x,
         haz       = "Hazard rate",
         cumhaz    = "Cumulative hazard rate",
         surv      = "Event free probability",
         cdf       = "Failure probability",
         loghaz    = "log(Hazard rate)",
         logcumhaz = "log(Cumulative hazard rate)",
         logsurv   = "log(Event free probability)",
         logcdf    = "log(Failure probability)",
         stop("Bug found: invalid input to 'type' argument."))
}

# Return appropriate y-axis limits for the prediction type
get_survpred_ylim <- function(x) {
  switch(x,
         surv = c(0,1),
         cdf  = c(0,1),
         NULL)
}

# Default plotting attributes
.PP_FILL <- "skyblue"
.PP_DARK <- "skyblue4"
.PP_VLINE_CLR <- "#222222"
.PP_YREP_CLR <- "#487575"
.PP_YREP_FILL <- "#222222"

# ----------------- internal -----------------------------------------------

handle_last_time <- function(object, last_time){
  if(is.null(last_time)){
    last_time <- lapply(seq_along(object$basehaz), function(b) NULL)
  } else {
    if(length(last_time) == 1){
      last_time <-  lapply(seq_along(object$basehaz), function(b) last_time)
    }
    if(all(seq_along(object$basehaz) != seq_along(last_time)))
      stop2("last_time has to be provided for each transition")
  }
  last_time
}


handle_times <- function(object, times){
  if(is.null(times)){
    times <- lapply(seq_along(object$basehaz), function(b) NULL)
  } else {
    if(length(times) == 1){
      times <-  lapply(seq_along(object$basehaz), function(b) times)
    }
    if(all(seq_along(object$basehaz) != seq_along(times)))
      stop2("times have to be provided for each transition")
  }
  times
}


handle_condition <- function(object, condition){
  if(is.null(condition)){
    condition <- lapply(seq_along(object$basehaz), function(b) NULL)
  } else {
    if(length(condition) == 1){
      condition <-  lapply(seq_along(object$basehaz), function(b) condition)
    }
    else if(all(seq_along(object$basehaz) != seq_along(condition)))
      stop2("condition should be provided for all transition")
  }
  condition
}

# Throw error if object isn't a idm object
#
# @param x The object to test.
handle_newdata <- function(object, newdata, call. = FALSE) {
  if (all(seq_along(object$basehaz) != seq_along(newdata)))
    stop("A new dataframe has to be provided for each transition.", call. = call.)
}



