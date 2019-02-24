#' Information criteria and cross-validation
#'
#' @description For models fit using MCMC, compute approximate leave-one-out
#'   cross-validation (LOO, LOOIC) or, less preferably, the Widely Applicable
#'   Information Criterion (WAIC) using the \pkg{\link[=loo-package]{loo}}
#'   package. Functions for \eqn{K}-fold cross-validation, model comparison,
#'   and model weighting/averaging are also provided. \strong{Note}:
#'   these functions are not guaranteed to work properly unless the \code{data}
#'   argument was specified when the model was fit. Also, as of \pkg{loo}
#'   version \code{2.0.0} the default number of cores is now only 1,  but we
#'   recommend using as many (or close to as many) cores as possible by setting
#'   the \code{cores} argument or using \code{options(mc.cores = VALUE)} to set
#'   it for an entire session.
#'
#' @aliases loo waic
#'
#' @export
#'
#' @param x For \code{loo}, \code{waic}, and \code{kfold} methods, a fitted
#'   model object returned by one of the rstanarm modeling functions. See
#'   \link{stanreg-objects}.
#'
#'   For \code{loo_model_weights}, \code{x} should be a "stanreg_list"
#'   object, which is a list of fitted model objects created by
#'   \code{\link{stanreg_list}}.
#'
#' @param ... For \code{compare_models}, \code{...} should contain two or more
#'   objects returned by the \code{loo}, \code{kfold}, or \code{waic} method
#'   (see the \strong{Examples} section, below).
#'
#'   For \code{loo_model_weights}, \code{...} should contain arguments
#'   (e.g. \code{method}) to pass to the default
#'   \code{\link[loo]{loo_model_weights}} method from the \pkg{loo} package.
#'
#' @param cores,save_psis Passed to \code{\link[loo]{loo}}.
#' @param k_threshold Threshold for flagging estimates of the Pareto shape
#'   parameters \eqn{k} estimated by \code{loo}. See the \emph{How to proceed
#'   when \code{loo} gives warnings} section, below, for details.
#'
#' @return The structure of the objects returned by \code{loo} and \code{waic}
#'   methods are documented in detail in the \strong{Value} section in
#'   \code{\link[loo]{loo}} and \code{\link[loo]{waic}} (from the \pkg{loo}
#'   package).
#'
#' @section Approximate LOO CV: The \code{loo} method for stanreg objects
#'   provides an interface to the \pkg{\link[=loo-package]{loo}} package for
#'   approximate leave-one-out cross-validation (LOO). The LOO Information
#'   Criterion (LOOIC) has the same purpose as the Akaike Information Criterion
#'   (AIC) that is used by frequentists. Both are intended to estimate the
#'   expected log predictive density (ELPD) for a new dataset. However, the AIC
#'   ignores priors and assumes that the posterior distribution is multivariate
#'   normal, whereas the functions from the \pkg{loo} package do not make this
#'   distributional assumption and integrate over uncertainty in the parameters.
#'   This only assumes that any one observation can be omitted without having a
#'   major effect on the posterior distribution, which can be judged using the
#'   diagnostic plot provided by the \code{\link[loo]{plot.loo}} method and the
#'   warnings provided by the \code{\link[loo]{print.loo}} method (see the
#'   \emph{How to Use the rstanarm Package} vignette for an example of this
#'   process).
#'
#'   \subsection{How to proceed when \code{loo} gives warnings (k_threshold)}{
#'   The \code{k_threshold} argument to the \code{loo} method for \pkg{rstanarm}
#'   models is provided as a possible remedy when the diagnostics reveal
#'   problems stemming from the posterior's sensitivity to particular
#'   observations. Warnings about Pareto \eqn{k} estimates indicate observations
#'   for which the approximation to LOO is problematic (this is described in
#'   detail in Vehtari, Gelman, and Gabry (2017) and the
#'   \pkg{\link[=loo-package]{loo}} package documentation). The
#'   \code{k_threshold} argument can be used to set the \eqn{k} value above
#'   which an observation is flagged. If \code{k_threshold} is not \code{NULL}
#'   and there are \eqn{J} observations with \eqn{k} estimates above
#'   \code{k_threshold} then when \code{loo} is called it will refit the
#'   original model \eqn{J} times, each time leaving out one of the \eqn{J}
#'   problematic observations. The pointwise contributions of these observations
#'   to the total ELPD are then computed directly and substituted for the
#'   previous estimates from these \eqn{J} observations that are stored in the
#'   object created by \code{loo}.
#'
#'   \strong{Note}: in the warning messages issued by \code{loo} about large
#'   Pareto \eqn{k} estimates we recommend setting \code{k_threshold} to at
#'   least \eqn{0.7}. There is a theoretical reason, explained in Vehtari,
#'   Gelman, and Gabry (2017), for setting the threshold to the stricter value
#'   of \eqn{0.5}, but in practice they find that errors in the LOO
#'   approximation start to increase non-negligibly when \eqn{k > 0.7}.
#'   }
#'
#' @seealso
#' \itemize{
#'   \item The new \href{http://mc-stan.org/loo/articles/}{\pkg{loo} package vignettes}
#'   and various \href{http://mc-stan.org/rstanarm/articles/}{\pkg{rstanarm} vignettes}
#'   for more examples using \code{loo} and related functions with \pkg{rstanarm} models.
#'   \item \code{\link[loo]{pareto-k-diagnostic}} in the \pkg{loo} package for
#'   more on Pareto \eqn{k} diagnostics.
#'   \item \code{\link{log_lik.stanreg}} to directly access the pointwise
#'   log-likelihood matrix.
#' }
#'
#' @examples
#' \donttest{
#' fit1 <- stan_glm(mpg ~ wt, data = mtcars)
#' fit2 <- stan_glm(mpg ~ wt + cyl, data = mtcars)
#'
#' # compare on LOOIC
#' # (for bigger models use as many cores as possible)
#' loo1 <- loo(fit1, cores = 2)
#' print(loo1)
#' loo2 <- loo(fit2, cores = 2)
#' print(loo2)
#'
#' # when comparing exactly two models, the reported 'elpd_diff'
#' # will be positive if the expected predictive accuracy for the
#' # second model is higher. the approximate standard error of the
#' # difference is also reported.
#' compare_models(loo1, loo2)
#' compare_models(loos = list(loo1, loo2)) # can also provide list
#'
#' # when comparing three or more models they are ordered by
#' # expected predictive accuracy. elpd_diff and se_diff are relative
#' # to the model with best elpd_loo (first row)
#' fit3 <- stan_glm(mpg ~ disp * as.factor(cyl), data = mtcars)
#' loo3 <- loo(fit3, cores = 2, k_threshold = 0.7)
#' compare_models(loo1, loo2, loo3)
#'
#' # setting detail=TRUE will also print model formulas
#' compare_models(loo1, loo2, loo3, detail=TRUE)
#'
#' # Computing model weights
#' model_list <- stanreg_list(fit1, fit2, fit3)
#' loo_model_weights(model_list, cores = 2) # can specify k_threshold=0.7 if necessary
#'
#' # if you have already computed loo then it's more efficient to pass a list
#' # of precomputed loo objects than a "stanreg_list", avoiding the need
#' # for loo_models weights to call loo() internally
#' loo_list <- list(fit1 = loo1, fit2 = loo2, fit3 = loo3) # names optional (affects printing)
#' loo_model_weights(loo_list)
#'
#' # 10-fold cross-validation
#' (kfold1 <- kfold(fit1, K = 10))
#' kfold2 <- kfold(fit2, K = 10)
#' compare_models(kfold1, kfold2, detail=TRUE)
#'
#' # Cross-validation stratifying by a grouping variable
#' # (note: might get some divergences warnings with this model but
#' # this is just intended as a quick example of how to code this)
#' library(loo)
#' fit4 <- stan_lmer(mpg ~ disp + (1|cyl), data = mtcars)
#' table(mtcars$cyl)
#' folds_cyl <- kfold_split_stratified(K = 3, x = mtcars$cyl)
#' table(cyl = mtcars$cyl, fold = folds_cyl)
#' kfold4 <- kfold(fit4, K = 3, folds = folds_cyl)
#' }
#'
#' @importFrom loo loo loo.function loo.matrix
#' @importFrom magrittr "%>%"
#'
loo.stanmstte <-
  function(x,
           ...,
           cores = getOption("mc.cores", 1),
           save_psis = FALSE,
           k_threshold = NULL) {
    if (!used.sampling(x))
      STOP_sampling_only("loo")
    if (model_has_weights(x))
      recommend_exact_loo(reason = "model has weights")

    user_threshold <- !is.null(k_threshold)
    if (user_threshold) {
      validate_k_threshold(k_threshold)
    } else {
      k_threshold <- 0.7
    }

    # chain_id to pass to loo::relative_eff
    chain_id <- chain_id_for_loo(x)

    has_quadrature <- any(x$has_quadrature)

    if (has_quadrature) {
      # ll <- log_lik.stanreg(x)
      # r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id, cores = cores)
      # loo_x <-
      #   suppressWarnings(loo.matrix(
      #     ll,
      #     r_eff = r_eff,
      #     cores = cores,
      #     save_psis = save_psis
      #   ))
      #  not implemented now.
    } else {
      args <- ll_args(x)
      llfun <- ll_fun(x)
      likfun <- function(data_i, draws) {
        exp(llfun(data_i, draws))
      }
      r_eff <- loo::relative_eff(
        # using function method
        x = likfun,
        chain_id = chain_id,
        data = args$data,
        draws = args$draws,
        cores = cores,
        ...
      )
      loo_x <- suppressWarnings(
        loo.function(
          llfun,
          data = args$data,
          draws = args$draws,
          r_eff = r_eff,
          ...,
          cores = cores,
          save_psis = save_psis
        )
      )
    }


    bad_obs <- loo::pareto_k_ids(loo_x, k_threshold)
    n_bad <- length(bad_obs)

    out <- structure(
      loo_x,
      name = deparse(substitute(x)),
      discrete = is_discrete(x),
      yhash = hash_y(x),
      formula = loo_model_formula(x)
    )

    if (!length(bad_obs)) {
      if (user_threshold) {
        message(
          "All pareto_k estimates below user-specified threshold of ",
          k_threshold,
          ". \nReturning loo object."
        )
      }
      return(out)
    }

    if (!user_threshold) {
      if (n_bad > 10) {
        recommend_kfold(n_bad)
      } else {
        recommend_reloo(n_bad)
      }
      return(out)
    }

    return(out)

    # reloo_out <- reloo(x, loo_x, obs = bad_obs)
    # structure(
    #   reloo_out,
    #   name = attr(out, "name"),
    #   discrete = attr(out, "discrete"),
    #   yhash = attr(out, "yhash"),
    #   formula = loo_model_formula(x)
    # )
    #
    # structure(
    #   reloo_out,
    #   name = attr(out, "name"),
    #   discrete = attr(out, "discrete"),
    #   yhash = attr(out, "yhash"),
    #   formula = loo_model_formula(x)
    # )

  }


# get arguments needed for ll_fun
# @param object stanmstte object
# @param newdata same as posterior predict
# @param offset vector of offsets (only required if model has offset term and
#   newdata is specified)
# @param m Integer specifying which submodel for stanmvreg objects
# @param reloo_or_kfold logical. TRUE if ll_args is for reloo or kfold
# @param ... For models without group-specific terms (i.e., not stan_[g]lmer),
#   if reloo_or_kfold is TRUE and 'newdata' is specified then ... is used to
#   pass 'newx' and 'stanmat' from reloo or kfold (bypassing pp_data). This is a
#   workaround in case there are issues with newdata containing factors with
#   only a single level. Or for stanmvreg objects, then ... can be used to pass
#   'stanmat', which may be a matrix with a reduced number of draws (potentially
#   just a single MCMC draw).
# @return a named list with elements data, draws, S (posterior sample size) and
#   N = number of observations
ll_args <- function(object, ...) UseMethod("ll_args")

#--- ll_args for stansurv models
ll_args.stanmstte <- function(object, newdata = NULL, ...) {

  validate_stanmstte_object(object)

  if (is.null(newdata)) {
    newdata <- get_model_data(object)
  }
  newdata <- lapply(newdata, function (nd) as.data.frame(nd))

  # response, ie. a Surv object
  form <- lapply(formula(object), function(f) as.formula(f) )
  y    <- lapply(seq_len(object$n_trans), function(n_t) eval(form[[n_t]][[2L]], newdata[[n_t]]) )

  # outcome, ie. time variables and status indicator
  t_beg   <- lapply(y, function(y_n) make_t(y_n, type = "beg") ) # entry time
  t_end   <- lapply(y, function(y_n) make_t(y_n,  type = "end") ) # exit  time
  t_upp   <- lapply(y, function(y_n) make_t(y_n,  type = "upp") ) # upper time for interval censoring
  status  <- lapply(y, function(y_n) make_d(y_n) )
  for(i in seq_along(status)){
    if (any(status[[i]] < 0 || status[[i]] > 3))
      stop2("Invalid status indicator in Surv object.")
  }

  # delayed entry indicator for each row of data
  delayed <- lapply(t_beg, function(t_beg_n) as.logical(!t_beg_n == 0) )

  # we reconstruct the design matrices even if no newdata, since it is
  # too much of a pain to store everything in the fitted model object
  # (e.g. w/ delayed entry, interval censoring, quadrature points, etc)
  pp <- pp_data(object, newdata, times = t_end)

  # returned object depends on quadrature not implemented
  if (any( object$has_quadrature) ){
    # pp_qpts_beg <- pp_data(object, newdata, times = t_beg, at_quadpoints = TRUE)
    # pp_qpts_end <- pp_data(object, newdata, times = t_end, at_quadpoints = TRUE)
    # pp_qpts_upp <- pp_data(object, newdata, times = t_upp, at_quadpoints = TRUE)
    # cpts <- c(pp$pts, pp_qpts_beg$pts, pp_qpts_end$pts, pp_qpts_upp$pts)
    # cwts <- c(pp$wts, pp_qpts_beg$wts, pp_qpts_end$wts, pp_qpts_upp$wts)
    # cids <- c(pp$ids, pp_qpts_beg$ids, pp_qpts_end$ids, pp_qpts_upp$ids)
    # x <- rbind(pp$x, pp_qpts_beg$x, pp_qpts_end$x, pp_qpts_upp$x)
    # s <- rbind(pp$s, pp_qpts_beg$s, pp_qpts_end$s, pp_qpts_upp$s)
    # x <- append_prefix_to_colnames(x, "x__")
    # s <- append_prefix_to_colnames(s, "s__")
    # status  <- c(status,  rep(NA, length(cids) - length(status)))
    # delayed <- c(delayed, rep(NA, length(cids) - length(delayed)))
    # data <- data.frame(cpts, cwts, cids, status, delayed)
    # data <- cbind(data, x, s)
  } else {
    x <- append_prefix_to_colnames(pp$x, "x__")
    cids <- lapply(newdata, function(d) get_id_var(d) )
    data <- dplyr::data_frame(
      cids = ulist(cids) ,
      t_beg = ulist(t_beg),
      t_end = ulist(t_end),
      t_upp = ulist(t_upp),
      status = ulist(status),
      delayed = ulist(delayed),
      type = uapply(seq_along(pp$ids), function(i) rep(i, length(pp$ids[[i]]))) )
    data <- dplyr::full_join(data, as.data.frame(x), by = c("cids" = "x__id_for_passing_to__"))
  }

  basehaz = lapply(seq_len(object$n_trans), function(i) get_basehaz(object, i) )

  # parameter draws
  draws                <- list()
  pars                 <- extract_pars(object)
  draws$basehaz        <- basehaz
  draws$n_trans        <- object$n_trans
  draws$transition_labels   <- object$transition_labels
  draws$aux            <- pars$aux
  draws$alpha          <- pars$alpha
  draws$beta           <- pars$beta
  draws$beta_tde       <- pars$beta_tde
  draws$has_quadrature <- pp$has_quadrature
  draws$qnodes         <- pp$qnodes

  out <- nlist(data, draws, S = NROW(draws$beta), N = n_distinct(cids))
  return(out)
}





#' Pointwise log-likelihood matrix
#'
#' For models fit using MCMC only, the \code{log_lik} method returns the
#' \eqn{S} by \eqn{N} pointwise log-likelihood matrix, where \eqn{S} is the size
#' of the posterior sample and \eqn{N} is the number of data points, or in the
#' case of the \code{stanmvreg} method (when called on \code{\link{stan_jm} and \code{\link{stan_mstte}}}
#' model objects) an \eqn{S} by \eqn{Npat} matrix where \eqn{Npat} is the number
#' of individuals.
#'
#' @aliases log_lik
#' @export
#'
#' @param newdata An optional data frame of new data (e.g. holdout data) to use
#'   when evaluating the log-likelihood. See the description of \code{newdata}
#'   for \code{\link{posterior_predict}}.
#' @param offset A vector of offsets. Only required if \code{newdata} is
#'   specified and an \code{offset} was specified when fitting the model.
#'
#' @return For the \code{stanreg} and \code{stanmvreg} methods an \eqn{S} by
#'   \eqn{N} matrix, where \eqn{S} is the size of the posterior sample and
#'   \eqn{N} is the number of data points. For the \code{stanjm} and \code{stanmstte} method
#'   an \eqn{S} by \eqn{Npat} matrix where \eqn{Npat} is the number of individuals.
#'
#'
#' @examples
#' \donttest{
#'  roaches$roach100 <- roaches$roach1 / 100
#'  fit <- stan_glm(
#'     y ~ roach100 + treatment + senior,
#'     offset = log(exposure2),
#'     data = roaches,
#'     family = poisson(link = "log"),
#'     prior = normal(0, 2.5),
#'     prior_intercept = normal(0, 10),
#'     iter = 500 # to speed up example
#'  )
#'  ll <- log_lik(fit)
#'  dim(ll)
#'  all.equal(ncol(ll), nobs(fit))
#'
#'  # using newdata argument
#'  nd <- roaches[1:2, ]
#'  nd$treatment[1:2] <- c(0, 1)
#'  ll2 <- log_lik(fit, newdata = nd, offset = c(0, 0))
#'  head(ll2)
#'  dim(ll2)
#'  all.equal(ncol(ll2), nrow(nd))
#' }
#'
log_lik.stanmstte <- function(object, newdata = NULL, ...) {
  # if (!used.sampling(object))
  #   STOP_sampling_only("Pointwise log-likelihood matrix")
  # validate_stanmstte_object(object)
  #
  # n_tra <- get_T(object)
  # if ("n_tra" %in% names(list(...)))
  #   stop("'n_tra' should not be specified for stan_mstte objects since the ",
  #        "log-likelihood is calculated for the full joint model.")
  # if(!is.null(newdata)) {
  #   newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
  #   if(is.data.frame(newdata)){
  #     if("trans" %!in% names(newdata) )
  #       stop2("New data-frame must be supplied with a column trans for transition.")
  #   } else if(is.list(newdata)){
  #     if(!identical(length(newdata), ntra) )
  #       stop2("Newdata must be supplied as a list with an element for each transition.")
  #   } else {
  #     stop2("Newdata must be supplied either as data-frame or list.")
  #   }
  # }
  #
  # if (!is.null(newdataLong)) {
  #   newdatas <- validate_newdatas(object, newdataLong, newdataEvent)
  #   newdataLong  <- newdatas[1:M]
  #   newdataEvent <- newdatas[["Event"]]
  # }
  # pars <- extract_pars(object) # full array of draws
  # data <- .pp_data_jm(object, newdataLong, newdataEvent)
  # calling_fun <- as.character(sys.call(-1))[1]
  # reloo_or_kfold <- calling_fun %in% c("kfold", "reloo")
  # val <- .ll_jm(object, data, pars, reloo_or_kfold = reloo_or_kfold, ...)
  # return(val)
  # not implemented
}

# Return a data frame for each submodel that:
# (1) only includes variables used in the model formula
# (2) only includes rows contained in the glmod/coxmod model frames
# (3) ensures that additional variables that are required
#     such as the ID variable or variables used in the
#     interaction-type association structures, are included.
#
# It is necessary to drop unneeded variables though so that
# errors are not encountered if the original data contained
# NA values for variables unrelated to the model formula.
# We generate a data frame here for in-sample predictions
# rather than using a model frame, since some quantities will
# need to be recalculated at quadrature points etc, for example
# in posterior_survfit.
#
# @param object A stanmstte object.
# @param m Integer specifying which submodel to get the
#   prediction data frame for (for stanmvreg or stanjm objects).
# @return A data frame or list of data frames with all the
#   (unevaluated) variables required for predictions.
get_model_data <- function(object, ...) UseMethod("get_model_data")

# ---- pp_data
get_model_data.stanmstte <- function(object, ...) {
  validate_stanmstte_object(object)
  terms <- terms(object)
  row_nms <- lapply(object$model_frame, function(mf) {
    row.names(mf)
  })
  lapply(seq_len(object$n_trans), function(n_t){
    newdata <- get_all_vars(terms[[n_t]], object$data[[n_t]])[row_nms[[n_t]], , drop = FALSE]
  })
}

get_model_data.stanmsjm <- function(object, ...) {
  validate_stanmsjm_object(object)
  M <- get_M(object)
  terms <- terms(object, fixed.only = FALSE)
  n_trans <- get_n_trans(object)
  transition_labels <- get_transition_name(object)
  
  extra_vars <- lapply(1:M, function(m) {
    # for each submodel loop over the four possible assoc  
    # interaction formulas and collect any variables used
    forms_m <- object$assoc["which_formulas",][[m]]
    uapply(forms_m, function(x) {
      if (length(x)) {
        rownames(attr(terms.formula(x), "factors")) 
      } else NULL
    })
  })
  # also ensure that id_var is in the event data
  for(n in seq_len(n_trans)){
    extra_vars[[paste0(transition_labels[[n]])]] <- object$id_var
  }
  n_extra_vars <- length(extra_vars) + n_trans - 1
  if (!identical(length(terms), length(extra_vars) ))
    stop2("Bug found: terms and extra_vars should be same length.")
  
  # add the extra variables to the terms formula for each submodel
  terms <- xapply(terms, extra_vars, FUN = function(x, y) {
    lhs <- x[[2L]]
    rhs <- deparse(x[[3L]], 500L)
    if (!is.null(y))
      rhs <- c(rhs, y)
    reformulate(rhs, response = lhs)
  })
  
  datas <- c(object$dataLong, list(object$dataEvent))

}

get_model_data <- function(object, m = NULL) {
  validate_stanmvreg_object(object)
  M <- get_M(object)
  terms <- terms(object, fixed.only = FALSE)
  
  # identify variables to add to the terms objects
  if (is.jm(object)) {
    extra_vars <- lapply(1:M, function(m) {
      # for each submodel loop over the four possible assoc  
      # interaction formulas and collect any variables used
      forms_m <- object$assoc["which_formulas",][[m]]
      uapply(forms_m, function(x) {
        if (length(x)) {
          rownames(attr(terms.formula(x), "factors")) 
        } else NULL
      })
    })
    # also ensure that id_var is in the event data
    extra_vars$Event <- object$id_var
    
    if (!identical(length(terms), length(extra_vars)))
      stop2("Bug found: terms and extra_vars should be same length.")
    
    # add the extra variables to the terms formula for each submodel
    terms <- xapply(terms, extra_vars, FUN = function(x, y) {
      lhs <- x[[2L]]
      rhs <- deparse(x[[3L]], 500L)
      if (!is.null(y))
        rhs <- c(rhs, y)
      reformulate(rhs, response = lhs)
    })
    
    datas <- c(object$dataLong, list(object$dataEvent))
  } else {
    datas <- object$data
  }
  
  # identify rows that were in the model frame
  row_nms <- lapply(model.frame(object), rownames)
  
  # drop rows and variables not required for predictions
  mfs <- xapply(w = terms, x = datas, y = row_nms,
                FUN = function(w, x, y) 
                  get_all_vars(w, x)[y, , drop = FALSE])
  
  mfs <- list_nms(mfs, M, stub = get_stub(object))
  if (is.null(m)) mfs else mfs[[m]]
}





pp_data <-
  function(object,
           newdata = NULL,
           re.form = NULL,
           offset = NULL,
           m = NULL,
           ...) {
    validate_stanmstte_object(object)
    if (is.stanmstte(object)) {
      return(.pp_data_mstte(object, newdata = newdata, ...))
    }
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
    qnodes <- 15L
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


# get log likelihood function for a particular model
# @param x stanreg object
# @return a function
ll_fun <- function(x, m = NULL) {
  validate_stanmstte_object(x)

  if (is.stansurv(x)) {
    return(.ll_surv_i)
  } else if (is.stanmstte(x) ){
    return(.ll_mstte_i)
  }

  fun <- paste0(".ll_", family(x, m = m)$family, "_i")
  get(fun, mode = "function")
}


#-----  Point likelihood for stanmstte

.ll_mstte_i <- function(data_i, draws) {

  if (any(draws$has_quadrature)) {

    # qnodes  <- draws$qnodes
    # status  <- data_i[1L, "status"]
    # delayed <- data_i[1L, "delayed"]
    #
    # # row indexing of quadrature points in data_i
    # idx_epts     <- 1
    # idx_qpts_beg <- 1 + (qnodes * 0) + (1:qnodes)
    # idx_qpts_end <- 1 + (qnodes * 1) + (1:qnodes)
    # idx_qpts_upp <- 1 + (qnodes * 2) + (1:qnodes)
    #
    # args <- list(times     = data_i$cpts,
    #              basehaz   = draws$basehaz,
    #              aux       = draws$aux,
    #              intercept = draws$alpha)
    #
    # eta  <- linear_predictor(draws$beta, .xdata_surv(data_i))
    # eta  <- eta + linear_predictor(draws$beta_tde, .sdata_surv(data_i))
    # lhaz <- eta + do.call(evaluate_log_basehaz, args)
    #
    # if (status == 1) {
    #   # uncensored
    #   lhaz_epts     <- lhaz[, idx_epts,     drop = FALSE]
    #   lhaz_qpts_end <- lhaz[, idx_qpts_end, drop = FALSE]
    #   lsurv <- -quadrature_sum(exp(lhaz_qpts_end),
    #                            qnodes = qnodes,
    #                            qwts   = data_i$cwts[idx_qpts_end])
    #   ll <- lhaz_epts + lsurv
    # } else if (status == 0) {
    #   # right censored
    #   lhaz_qpts_end <- lhaz[, idx_qpts_end, drop = FALSE]
    #   lsurv <- -quadrature_sum(exp(lhaz_qpts_end),
    #                            qnodes = qnodes,
    #                            qwts   = data_i$cwts[idx_qpts_end])
    #   ll <- lsurv
    # } else if (status == 2) {
    #   # left censored
    #   lhaz_qpts_end <- lhaz[, idx_qpts_end, drop = FALSE]
    #   lsurv <- -quadrature_sum(exp(lhaz_qpts_end),
    #                            qnodes = qnodes,
    #                            qwts   = data_i$cwts[idx_qpts_end])
    #   ll <- log(1 - exp(lsurv)) # = log CDF
    # } else if (status == 3) {
    #   # interval censored
    #   lhaz_qpts_end <- lhaz[, idx_qpts_end, drop = FALSE]
    #   lsurv_lower <- -quadrature_sum(exp(lhaz_qpts_end),
    #                                  qnodes = qnodes,
    #                                  qwts   = data_i$cwts[idx_qpts_end])
    #   lhaz_qpts_upp <- lhaz[, idx_qpts_upp, drop = FALSE]
    #   lsurv_upper <- -quadrature_sum(exp(lhaz_qpts_upp),
    #                                  qnodes = qnodes,
    #                                  qwts   = data_i$cwts[idx_qpts_upp])
    #   ll <- log(exp(lsurv_lower) - exp(lsurv_upper))
    # }
    # if (delayed) {
    #   # delayed entry
    #   lhaz_qpts_beg <- lhaz[, idx_qpts_beg, drop = FALSE]
    #   lsurv_beg <- -quadrature_sum(exp(lhaz_qpts_beg),
    #                                qnodes = qnodes,
    #                                qwts   = data_i$cwts[idx_qpts_beg])
    #   ll <- ll - lsurv_beg
    # }

  } else { # no quadrature

    status  <- data_i$status
    delayed <- data_i$delayed

    nms_aux <- grep(get_transition_name(draws, data_i$type), colnames(draws$aux), value = TRUE, fixed = TRUE)
    nms_alpha <- grep(get_transition_name(draws, data_i$type), colnames(draws$alpha), value = TRUE, fixed = TRUE)

    args <- list(basehaz   = draws$basehaz[[data_i$type]],
                 aux       = draws$aux[, nms_aux, drop = FALSE],
                 intercept =  draws$alpha[, nms_alpha, drop = FALSE])

    eta  <- linear_predictor(draws$beta, .xdata_surv(data_i))

    if (status == 1) {
      # uncensored
      args$times <- data_i$t_end

      lhaz  <- do.call(evaluate_log_basehaz,  args) + eta
      lsurv <- do.call(evaluate_log_basesurv, args) * exp(eta)
      ll <- lhaz + lsurv
    } else if (status == 0) {
      # right censored
      args$times <- data_i$t_end
      lsurv <- do.call(evaluate_log_basesurv, args) * exp(eta)
      ll <- lsurv
    } else if (status == 2) {
      # left censored
      args$times <- data_i$t_end
      lsurv <- do.call(evaluate_log_basesurv, args) * exp(eta)
      ll <- log(1 - exp(lsurv)) # = log CDF
    } else if (status == 3) {
      # interval censored
      args$times  <- data_i$t_end
      lsurv_lower <- do.call(evaluate_log_basesurv, args) * exp(eta)
      args$times  <- data_i$t_upp
      lsurv_upper <- do.call(evaluate_log_basesurv, args) * exp(eta)
      ll <- log(exp(lsurv_lower) - exp(lsurv_upper))
    }
    if (delayed) {
      # delayed entry
      args$times <- data_i$t_beg
      lsurv_beg <- do.call(evaluate_log_basesurv, args) * exp(eta)
      ll <- ll - lsurv_beg
    }

  }
  return(ll)
}

#------------- internal
#------------------------------
# Below are code chunks taken from the 'rstanarm' R package, obtained
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University


recommend_exact_loo <- function(reason) {
  stop(
    "'loo' is not supported if ", reason, ". ",
    "If refitting the model 'nobs(x)' times is feasible, ",
    "we recommend calling 'kfold' with K equal to the ",
    "total number of observations in the data to perform exact LOO-CV.\n",
    call. = FALSE
  )
}

# chain_id to pass to loo::relative_eff
chain_id_for_loo <- function(object) {
  dims <- dim(object$stanfit)[1:2]
  n_iter <- dims[1]
  n_chain <- dims[2]
  rep(1:n_chain, each = n_iter)
}


# for each transition
.xdata_surv <- function(data) {
  nms <- colnames(data)
  sel <- grep("^x__", nms)
  data[, sel]
}

# Evaluate the log baseline hazard at the specified times given the
# vector or matrix of MCMC draws for the baseline hazard parameters
#
# @param times A vector of times.
# @param basehaz A list with info about the baseline hazard.
# @param aux,intercept A vector or matrix of parameter estimates (MCMC draws).
# @param x Predictor matrix.
# @param s Predictor matrix for time-dependent effects.
# @return A vector or matrix, depending on the input type of aux.
evaluate_log_basehaz <- function(times, basehaz, aux, intercept = NULL) {
  switch(get_basehaz_name(basehaz),
         "exp"       = log_basehaz_exponential(times, log_scale = intercept),
         "weibull"   = log_basehaz_weibull (times, shape = aux, log_scale = intercept),
         "gompertz"  = log_basehaz_gompertz(times, scale = aux, log_shape = intercept),
         "ms"        = log_basehaz_ms(times, coefs = aux, basis = basehaz$basis),
         "bs"        = log_basehaz_bs(times, coefs = aux, basis = basehaz$basis),
         "piecewise" = log_basehaz_pw(times, coefs = aux, knots = basehaz$knots),
         stop2("Bug found: unknown type of baseline hazard."))
}

log_basehaz_exponential <- function(x, log_scale) {
  linear_predictor(log_scale, rep(1, length(x)))
}
log_basehaz_weibull  <- function(x, shape, log_scale) {
  as.vector(log_scale + log(shape)) + linear_predictor(shape - 1, log(x))
}
log_basehaz_gompertz <- function(x, scale, log_shape) {
  as.vector(log_shape) + linear_predictor(scale, x)
}
log_basehaz_ms <- function(x, coefs, basis) {
  log(linear_predictor(coefs, basis_matrix(x, basis = basis)))
}
log_basehaz_bs <- function(x, coefs, basis) {
  linear_predictor(coefs, basis_matrix(x, basis = basis))
}
log_basehaz_pw <- function(x, coefs, knots) {
  linear_predictor(coefs, dummy_matrix(x, knots = knots))
}

# Evaluate the log baseline hazard at the specified times
# given the vector or matrix of MCMC draws for the baseline
# hazard coeffients / parameters
#
# @param log_haz A vector containing the log hazard for each
#   individual, evaluated at each of the quadrature points. The
#   vector should be ordered such that the first N elements contain
#   the log_haz evaluated for each individual at quadrature point 1,
#   then the next N elements are the log_haz evaluated for each
#   individual at quadrature point 2, and so on.
# @param qnodes Integer specifying the number of quadrature nodes
#   at which the log hazard was evaluated for each individual.
# @param qwts A vector of unstandardised GK quadrature weights.
# @return A vector or matrix of log survival probabilities.
evaluate_log_survival <- function(log_haz, qnodes, qwts) {
  UseMethod("evaluate_log_survival")
}

evaluate_log_survival.default <- function(log_haz, qnodes, qwts) {
  # convert log hazard to hazard
  haz <- exp(log_haz)
  # apply GK quadrature weights
  weighted_haz <- qwts * haz
  # sum quadrature points for each individual to get cum_haz
  splitting_vec <- rep(1:qnodes, each = length(haz) / qnodes)
  cum_haz <- Reduce('+', split(weighted_haz, splitting_vec))
  # return: -cum_haz == log survival probability
  -cum_haz
}

evaluate_log_survival.matrix <- function(log_haz, qnodes, qwts) {
  # convert log hazard to hazard
  haz <- exp(log_haz)
  # apply GK quadrature weights
  weighted_haz <- sweep(haz, 2L, qwts, `*`)
  # sum quadrature points for each individual to get cum_haz
  cum_haz <- Reduce('+', array2list(weighted_haz, nsplits = qnodes))
  # return: -cum_haz == log survival probability
  -cum_haz
}

# Evaluate the log baseline hazard at the specified times
# given the vector or matrix of MCMC draws for the baseline
# hazard coeffients / parameters
#
# @param times A vector of times.
# @param basehaz A list with info about the baseline hazard.
# @param coefs A vector or matrix of parameter estimates (MCMC draws).
# @return A vector or matrix, depending on the input type of coefs.
evaluate_log_basehaz2 <- function(times, basehaz, coefs) {
  type <- basehaz$type_name
  if (type == "weibull") {
    X  <- log(times) # log times
    B1 <- log(coefs) # log shape
    B2 <- coefs - 1  # shape - 1
    log_basehaz <- as.vector(B1) + linear_predictor(B2,X)
  } else if (type == "bs") {
    X <- predict(basehaz$bs_basis, times) # b-spline basis
    B <- coefs                            # b-spline coefs
    log_basehaz <- linear_predictor(B,X)
  } else {
    stop2("Not yet implemented for basehaz = ", type)
  }
  log_basehaz
}



evaluate_log_haz <- function(times, basehaz, betas, betas_tde, aux,
                             intercept = NULL, x, s = NULL) {
  eta <- linear_predictor(betas, x)
  if ((!is.null(s)) && ncol(s))
    eta <- eta + linear_predictor(betas_tde, s)
  args <- nlist(times, basehaz, aux, intercept)
  do.call(evaluate_log_basehaz, args) + eta
}

evaluate_basehaz <- function(times, basehaz, aux, intercept = NULL) {
  exp(evaluate_log_basehaz(times = times, basehaz = basehaz,
                           aux = aux, intercept = intercept))
}

#-------------

# Evaluate the log baseline survival at the specified times given the
# vector or matrix of MCMC draws for the baseline hazard parameters
#
# @param times A vector of times.
# @param basehaz A list with info about the baseline hazard.
# @param aux,intercept A vector or matrix of parameter estimates (MCMC draws).
# @return A vector or matrix, depending on the input type of aux.
evaluate_log_basesurv <- function(times, basehaz, aux, intercept = NULL) {
  switch(get_basehaz_name(basehaz),
         "exp"       = log_basesurv_exponential(times, log_scale = intercept),
         "weibull"   = log_basesurv_weibull (times, shape = aux, log_scale = intercept),
         "gompertz"  = log_basesurv_gompertz(times, scale = aux, log_shape = intercept),
         "ms"        = log_basesurv_ms(times, coefs = aux, basis = basehaz$basis),
         stop2("Bug found: unknown type of baseline hazard."))
}

log_basesurv_exponential <- function(x, log_scale) {
  -linear_predictor(exp(log_scale), x)
}
log_basesurv_weibull  <- function(x, shape, log_scale) {
  -exp(as.vector(log_scale) + linear_predictor(shape, log(x)))
}
log_basesurv_gompertz <- function(x, scale, log_shape) {
  -(as.vector(exp(log_shape) / scale)) * (exp(linear_predictor(scale, x)) - 1)
}
log_basesurv_ms <- function(x, coefs, basis) {
  -linear_predictor(coefs, basis_matrix(x, basis = basis, integrate = TRUE))
}

evaluate_log_surv <- function(times, basehaz, betas, aux, intercept = NULL, x, ...) {
  eta  <- linear_predictor(betas, x)
  args <- nlist(times, basehaz, aux,  intercept)
  do.call(evaluate_log_basesurv, args) * exp(eta)
}

#---------------

quadrature_sum <- function(x, qnodes, qwts) {
  UseMethod("quadrature_sum")
}

quadrature_sum.default <- function(x, qnodes, qwts) {
  weighted_x <- qwts * x                                 # apply quadrature weights
  splitted_x <- split_vector(x, n_segments = qnodes)     # split at each quad node
  Reduce('+', splitted_x)                                # sum over the quad nodes
}

quadrature_sum.matrix <- function(x, qnodes, qwts) {
  weighted_x <- sweep_multiply(x, qwts, margin = 2L)     # apply quadrature weights
  splitted_x <- array2list(weighted_x, nsplits = qnodes) # split at each quad node
  Reduce('+', splitted_x)                                # sum over the quad nodes
}
