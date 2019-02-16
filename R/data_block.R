#------------------------------
# Below are code chunks taken from the 'rstanarm' R package, obtained
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# Copyright (C) 2016, 2017 Sam Brilleman

#' @importFrom survival Surv
#' @export
survival::Surv


# Check input argument is a valid type, and return as a list
#
# @param arg The user input to the argument
# @param type A character vector of valid classes
# @param validate_length The required length of the returned list
# @return A list
validate_arg <- function(arg, type, validate_length = NULL) {
  nm <- deparse(substitute(arg))

  if (inherits(arg, type)) {
    # input type is valid, so return as a list
    arg <- list(arg)
  }
  else if (is(arg, "list")) {
    # input type is a list, check each element
    check <- sapply(arg, function(x) inherits(x, type))
    if (!all(check))
      STOP_arg(nm, type)
  }
  else {
    # input type is not valid
    STOP_arg(nm, type)
  }

  if (!is.null(validate_length)) {
    # return list of the specified length
    if (length(arg) == 1L)
      arg <- rep(arg, times = validate_length)
    if (!length(arg) == validate_length)
      stop2(nm, " is a list of the incorrect length.")
  }

  if ("data.frame" %in% type)
    arg <- lapply(arg, as.data.frame)
  if ("family" %in% type)
    arg <- lapply(arg, validate_family)

  arg
}


# Center a matrix x and return extra stuff
#
# @param x A design matrix
# @param sparse A flag indicating whether x is to be treated as sparse
center_x <- function(x, sparse) {
  x <- as.matrix(x)
  has_intercept <- if (ncol(x) == 0)
    FALSE else grepl("(Intercept", colnames(x)[1L], fixed = TRUE)

  xtemp <- if (has_intercept) x[, -1L, drop=FALSE] else x
  if (has_intercept && !sparse) {
    xbar <- colMeans(xtemp)
    xtemp <- sweep(xtemp, 2, xbar, FUN = "-")
  }
  else xbar <- rep(0, ncol(xtemp))

  sel <- apply(xtemp, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  if (any(sel)) {
    # drop any column of x with < 2 unique values (empty interaction levels)
    # exception is column of 1s isn't dropped
    warning("Dropped empty interaction levels: ",
            paste(colnames(xtemp)[sel], collapse = ", "))
    xtemp <- xtemp[, !sel, drop = FALSE]
    xbar <- xbar[!sel]
  }

  return(nlist(xtemp, xbar, has_intercept))
}


# Function to return the range or SD of the predictors, used for scaling the priors
# This is taken from an anonymous function in stan_glm.fit
#
# @param x A vector
get_scale_value <- function(x) {
  num.categories <- dplyr::n_distinct(x)
  x.scale <- 1
  if (num.categories == 2) {
    x.scale <- diff(range(x))
  } else if (num.categories > 2) {
    x.scale <- sd(x)
  }
  return(x.scale)
}

# Deal with priors
#
# @param prior A list
# @param nvars An integer indicating the number of variables
# @param default_scale Default value to use to scale if not specified by user
# @param link String naming the link function.
# @param ok_dists A list of admissible distributions.
handle_glm_prior <- function(prior, nvars, default_scale, link,
                             ok_dists = nlist("normal", student_t = "t",
                                              "cauchy", "hs", "hs_plus",
                                              "laplace", "lasso", "product_normal")) {
  if (!length(prior))
    return(list(prior_dist = 0L, prior_mean = as.array(rep(0, nvars)),
                prior_scale = as.array(rep(1, nvars)),
                prior_df = as.array(rep(1, nvars)), prior_dist_name = NA,
                global_prior_scale = 0, global_prior_df = 0,
                slab_df = 0, slab_scale = 0,
                prior_autoscale = FALSE))

  if (!is.list(prior))
    stop(sQuote(deparse(substitute(prior))), " should be a named list")

  prior_dist_name <- prior$dist
  prior_scale <- prior$scale
  prior_mean <- prior$location
  prior_df <- prior$df
  prior_mean[is.na(prior_mean)] <- 0
  prior_df[is.na(prior_df)] <- 1
  global_prior_scale <- 0
  global_prior_df <- 0
  slab_df <- 0
  slab_scale <- 0
  if (!prior_dist_name %in% unlist(ok_dists)) {
    stop("The prior distribution should be one of ",
         paste(names(ok_dists), collapse = ", "))
  } else if (prior_dist_name %in%
             c("normal", "t", "cauchy", "laplace", "lasso", "product_normal")) {
    if (prior_dist_name == "normal") prior_dist <- 1L
    else if (prior_dist_name == "t") prior_dist <- 2L
    else if (prior_dist_name == "laplace") prior_dist <- 5L
    else if (prior_dist_name == "lasso") prior_dist <- 6L
    else if (prior_dist_name == "product_normal") prior_dist <- 7L
    prior_scale <- set_prior_scale(prior_scale, default = default_scale,
                                   link = link)
  } else if (prior_dist_name %in% c("hs", "hs_plus")) {
    prior_dist <- ifelse(prior_dist_name == "hs", 3L, 4L)
    global_prior_scale <- prior$global_scale
    global_prior_df <- prior$global_df
    slab_df <- prior$slab_df
    slab_scale <- prior$slab_scale
  } else if (prior_dist_name %in% "exponential") {
    prior_dist <- 3L # only used for scale parameters so 3 not a conflict with 3 for hs
  }

  prior_df <- maybe_broadcast(prior_df, nvars)
  prior_df <- as.array(pmin(.Machine$double.xmax, prior_df))
  prior_mean <- maybe_broadcast(prior_mean, nvars)
  prior_mean <- as.array(prior_mean)
  prior_scale <- maybe_broadcast(prior_scale, nvars)

  nlist(prior_dist,
        prior_mean,
        prior_scale,
        prior_df,
        prior_dist_name,
        global_prior_scale,
        global_prior_df,
        slab_df,
        slab_scale,
        prior_autoscale = isTRUE(prior$autoscale))
}


# Check and set scale parameters for priors
#
# @param scale Value of scale parameter (can be NULL).
# @param default Default value to use if \code{scale} is NULL.
# @param link String naming the link function or NULL.
# @return If a probit link is being used, \code{scale} (or \code{default} if
#   \code{scale} is NULL) is scaled by \code{dnorm(0) / dlogis(0)}. Otherwise
#   either \code{scale} or \code{default} is returned.
set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link) || is.null(link))
  if (is.null(scale))
    scale <- default
  if (isTRUE(link == "probit"))
    scale <- scale * dnorm(0) / dlogis(0)

  return(scale)
}

# Autoscaling of priors
#
# @param prior_stuff A named list returned by a call to handle_glm_prior
# @param response A vector containing the response variable, only required if
#   the priors are to be scaled by the standard deviation of the response (for
#   gaussian reponse variables only)
# @param predictors The predictor matrix, only required if the priors are to be
#   scaled by the range/sd of the predictors
# @param family A family object
# @param QR A logical specifying whether QR decomposition is used for the
#   predictor matrix
# @param min_prior_scale The minimum allowed for prior scales
# @param assoc A two dimensional array with information about desired association
#   structure for the joint model (returned by a call to validate_assoc). Cannot
#   be NULL if autoscaling priors for the association parameters.
# @param ... Other arguments passed to make_assoc_terms. If autoscaling priors
#   for the association parameters then this should include 'parts' which
#   is a list containing the design matrices for the longitudinal submodel
#   evaluated at the quadrature points, as well as 'beta' and 'b' which are
#   the parameter values to use when constructing the linear predictor(s) in
#   make_assoc_terms.
# @return A named list with the same structure as returned by handle_glm_prior
autoscale_prior <- function(prior_stuff, response = NULL, predictors = NULL,
                            family = NULL, QR = FALSE, min_prior_scale = 1e-12,
                            assoc = NULL, ...) {
  ps <- prior_stuff

  if (!is.null(response) && is.gaussian(family)) {
    # use response variable for scaling priors
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      ss <- sd(response)
      ps$prior_scale <- ss * ps$prior_scale
    }
  }

  if (!is.null(predictors) && !QR) {
    # use predictors for scaling priors
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      ps$prior_scale <-
        pmax(min_prior_scale,
             ps$prior_scale / apply(predictors, 2L, get_scale_value))
    }
  }

  if (!is.null(assoc)) {
    # Evaluate mean and SD of each of the association terms that will go into
    # the linear predictor for the event submodel (as implicit "covariates").
    # (NB the approximate association terms are calculated using coefs
    # from the separate longitudinal submodels estimated using glmer).
    # The mean will be used for centering each association term.
    # The SD will be used for autoscaling the prior for each association parameter.
    if (is.null(family))
      stop("'family' cannot be NULL when autoscaling association parameters.")
    assoc_terms <- make_assoc_terms(family = family, assoc = assoc, ...)
    ps$a_xbar <- as.array(apply(assoc_terms, 2L, mean))
    if (ps$prior_dist > 0L && ps$prior_autoscale) {
      a_beta_scale <- apply(assoc_terms, 2L, get_scale_value)
      ps$prior_scale <- pmax(min_prior_scale, ps$prior_scale / a_beta_scale)
    }
  }

  ps$prior_scale <- as.array(pmin(.Machine$double.xmax, ps$prior_scale))
  ps
}

# Return the default scale parameter for 'prior_aux'.
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A scalar.
get_default_aux_scale <- function(basehaz) {
  nm <- get_basehaz_name(basehaz)
  if (nm %in% c("weibull", "gompertz")) 2 else 20
}

# Create "prior.info" attribute for stan_{mvmer,jm}; needed for prior_summary()
#
# @param user_* The user's priors. These should be passed in after broadcasting
#   the df/location/scale arguments if necessary.
# @param y_has_intercept Vector of T/F, does each long submodel have an intercept?
# @param y_has_predictors Vector of T/F, does each long submodel have predictors?
# @param e_has_intercept T/F, does event submodel have an intercept?
# @param e_has_predictors T/F, does event submodel have predictors?
# @param has_assoc Logical specifying whether the model has an association
#   structure. Can be NULL if the prior summary is not for a joint model.
# @param adjusted_prior_*_scale Adjusted scales computed if using autoscaled priors
# @param family A list of family objects.
# @param basehaz A list with information about the baseline hazard.
# @param stub_for_names Character string with the text stub to use in the
#   names identifying the glmer or longitudinal submodels.
# @return A named list with components 'prior*', 'prior*_intercept',
#   'prior_covariance' and 'prior*_aux' each of which itself is a list
#   containing the needed values for prior_summary.
summarize_jm_prior <-
  function(user_priorLong = NULL,
           user_priorLong_intercept = NULL,
           user_priorLong_aux = NULL,
           user_priorEvent = NULL,
           user_priorEvent_intercept = NULL,
           user_priorEvent_aux = NULL,
           user_priorEvent_assoc = NULL,
           user_prior_covariance = NULL,
           b_user_prior_stuff = NULL,
           b_prior_stuff = NULL,
           y_has_intercept = NULL,
           e_has_intercept = NULL,
           y_has_predictors = NULL,
           e_has_predictors = NULL,
           has_assoc = NULL,
           adjusted_priorLong_scale = NULL,
           adjusted_priorLong_intercept_scale = NULL,
           adjusted_priorLong_aux_scale = NULL,
           adjusted_priorEvent_scale = NULL,
           adjusted_priorEvent_intercept_scale = NULL,
           adjusted_priorEvent_aux_scale = NULL,
           adjusted_priorEvent_assoc_scale = NULL,
           family = NULL,
           basehaz = NULL,
           stub_for_names = "Long") {
    if (!is.null(family) && !is(family, "list"))
      stop("'family' should be a list of family objects, one for each submodel.")
    if (!is.null(has_assoc) && !is.logical(has_assoc) && (length(has_assoc) == 1L))
      stop("'has_assoc' should be a logical vector of length 1.")
    M <- length(family)

    prior_list <- list()

    if (!is.null(user_priorLong)) {
      rescaled_coefLong <- mapply(check_if_rescaled, user_priorLong,
                                  y_has_predictors, adjusted_priorLong_scale)
      rescaled_intLong  <- mapply(check_if_rescaled, user_priorLong_intercept,
                                  y_has_intercept, adjusted_priorLong_intercept_scale)
      rescaled_auxLong  <- mapply(check_if_rescaled, user_priorLong_aux,
                                  TRUE, adjusted_priorLong_aux_scale)
      for (m in 1:M) {
        user_priorLong[[m]] <-
          rename_t_and_cauchy(user_priorLong[[m]], y_has_predictors[m])
        user_priorLong_intercept[[m]] <-
          rename_t_and_cauchy(user_priorLong_intercept[[m]], y_has_intercept[m])
        user_priorLong_aux[[m]] <-
          rename_t_and_cauchy(user_priorLong_aux[[m]], TRUE)
      }
      prior_list$priorLong <- list_nms(lapply(1:M, function(m) {
        if (!y_has_predictors[m]) NULL else with(user_priorLong[[m]], list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coefLong[m])
            adjusted_priorLong_scale[[m]] else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        ))
      }), M, stub = stub_for_names)
      prior_list$priorLong_intercept <- list_nms(lapply(1:M, function(m) {
        if (!y_has_intercept[m]) NULL else with(user_priorLong_intercept[[m]], list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_intLong[m])
            adjusted_priorLong_intercept_scale[[m]] else NULL,
          df = if (prior_dist_name %in% "student_t")
            prior_df else NULL
        ))
      }), M, stub = stub_for_names)
      aux_name <- lapply(family, .rename_aux)
      prior_list$priorLong_aux <- list_nms(lapply(1:M, function(m) {
        if (is.na(aux_name[[m]])) NULL else with(user_priorLong_aux[[m]], list(
          dist = prior_dist_name,
          location = if (!is.na(prior_dist_name) &&
                         prior_dist_name != "exponential")
            prior_mean else NULL,
          scale = if (!is.na(prior_dist_name) &&
                      prior_dist_name != "exponential")
            prior_scale else NULL,
          adjusted_scale = if (rescaled_auxLong[m])
            adjusted_priorLong_aux_scale[[m]] else NULL,
          df = if (!is.na(prior_dist_name) &&
                   prior_dist_name %in% "student_t")
            prior_df else NULL,
          rate = if (!is.na(prior_dist_name) &&
                     prior_dist_name %in% "exponential")
            1 / prior_scale else NULL,
          aux_name = aux_name[[m]]
        ))
      }), M, stub = stub_for_names)
    }

    if (!is.null(user_priorEvent)) {
      rescaled_coefEvent <- check_if_rescaled(user_priorEvent, e_has_predictors,
                                              adjusted_priorEvent_scale)
      rescaled_intEvent  <- check_if_rescaled(user_priorEvent_intercept, e_has_intercept,
                                              adjusted_priorEvent_intercept_scale)
      rescaled_auxEvent  <- check_if_rescaled(user_priorEvent_aux, TRUE,
                                              adjusted_priorEvent_aux_scale)
      user_priorEvent <-
        rename_t_and_cauchy(user_priorEvent, e_has_predictors)
      user_priorEvent_intercept <-
        rename_t_and_cauchy(user_priorEvent_intercept, e_has_intercept)
      user_priorEvent_aux <-
        rename_t_and_cauchy(user_priorEvent_aux, TRUE)
      prior_list$priorEvent <-
        if (!e_has_predictors) NULL else with(user_priorEvent, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coefEvent)
            adjusted_priorEvent_scale else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        ))
      prior_list$priorEvent_intercept <-
        if (!e_has_intercept) NULL else with(user_priorEvent_intercept, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_intEvent)
            adjusted_priorEvent_intercept_scale else NULL,
          df = if (prior_dist_name %in% "student_t")
            prior_df else NULL
        ))
      e_aux_name <- .rename_e_aux(basehaz)
      prior_list$priorEvent_aux <-
        with(user_priorEvent_aux, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_auxEvent)
            adjusted_priorEvent_aux_scale else NULL,
          df = if (!is.na(prior_dist_name) &&
                   prior_dist_name %in% "student_t")
            prior_df else NULL,
          aux_name = e_aux_name
        ))
    }

    if (!is.null(user_priorEvent_assoc)) {
      rescaled_coefAssoc <- check_if_rescaled(user_priorEvent_assoc, has_assoc,
                                              adjusted_priorEvent_assoc_scale)
      user_priorEvent_assoc <- rename_t_and_cauchy(user_priorEvent_assoc, has_assoc)
      prior_list$priorEvent_assoc <-
        if (!has_assoc) NULL else with(user_priorEvent_assoc, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coefAssoc)
            adjusted_priorEvent_assoc_scale else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        ))
    }

    if (length(user_prior_covariance)) {
      if (user_prior_covariance$dist == "decov") {
        prior_list$prior_covariance <- user_prior_covariance
      } else if (user_prior_covariance$dist == "lkj") {
        # lkj prior for correlation matrix
        prior_list$prior_covariance <- user_prior_covariance
        # half-student_t prior on SD for each ranef (possibly autoscaled)
        prior_list$prior_covariance$df <- b_user_prior_stuff$prior_df
        prior_list$prior_covariance$scale <- b_user_prior_stuff$prior_scale
        adj_scales <- uapply(b_prior_stuff, FUN = uapply, '[[', "prior_scale")
        if (!all(b_user_prior_stuff$prior_scale == adj_scales)) {
          prior_list$prior_covariance$adjusted_scale <- adj_scales
        } else {
          prior_list$prior_covariance$adjusted_scale <- NULL
        }
      } else {
        prior_list$prior_covariance <- NULL
      }
    }

    if (!stub_for_names == "Long") {
      nms <- names(prior_list)
      new_nms <- gsub("Long", "", nms)
      names(prior_list) <- new_nms
    }

    return(prior_list)
  }


# Check the family and link function are supported by stan_{mvmer,jm}
#
# @param family A family object
# @param supported_families A character vector of supported family names
# @return A family object
validate_famlink <- function(family, supported_families) {
  famname <- family$family
  fam <- which(supported_families == famname)
  if (!length(fam))
    stop2("'family' must be one of ", paste(supported_families, collapse = ", "))
  supported_links <- supported_glm_links(famname)
  link <- which(supported_links == family$link)
  if (!length(link))
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  return(family)
}

# @param famname string naming the family
# @return character vector of supported link functions for the family
supported_glm_links <- function(famname) {
  switch(
    famname,
    binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
    gaussian = c("identity", "log", "inverse"),
    Gamma = c("identity", "log", "inverse"),
    inverse.gaussian = c("identity", "log", "inverse", "1/mu^2"),
    "neg_binomial_2" = , # intentional
    poisson = c("log", "identity", "sqrt"),
    "Beta regression" = c("logit", "probit", "cloglog", "cauchit"),
    stop("unsupported family")
  )
}


# Check if the user input a list of priors for the longitudinal
# submodel, and if not, then return the appropriate list
#
# @param prior The user input to the prior argument in the stan_mvmer
#   or stan_jm call
# @param M An integer specifying the number of longitudinal submodels
broadcast_prior <- function(prior, M) {
  if (is.null(prior)) {
    return(rep(list(NULL), M))
  }
  else if ("dist" %in% names(prior)) {
    return(rep(list(prior), M))
  }
  else if (is.list(prior) && length(prior) == M) {
    return(prior)
  }
  else {
    nm <- deparse(substitute(priorarg))
    stop2(nm, " appears to provide prior information separately for the ",
          "different submodels, but the list is of the incorrect length.")
  }
}

# Construct a list with information about the longitudinal submodel
#
# @param formula The model formula for the glmer submodel.
# @param data The data for the glmer submodel.
# @param family The family object for the glmer submodel.
# @return A named list with the following elements:
#   y: named list with the reponse vector and related info.
#   x: named list with the fe design matrix and related info.
#   z: named list with the re design matrices and related info.
#   terms: the model.frame terms object with bars "|" replaced by "+".
#   model_frame: The model frame with all variables used in the
#     model formula.
#   formula: The model formula.
#   reTrms: returned by lme4::glFormula$reTrms.
#   family: the (modified) family object for the glmer submodel.
#   intercept_type: named list with info about the type of
#     intercept required for the glmer submodel.
#   has_aux: logical specifying whether the glmer submodel
#     requires an auxiliary parameter.
handle_y_mod <- function(formula, data, family, stub) {
  mf <- stats::model.frame(lme4::subbars(formula), data)
  if (!length(formula) == 3L)
    stop2("An outcome variable must be specified.")

  # lme4 parts
  lme4_parts <- lme4::glFormula(formula, data)
  reTrms <- lme4_parts$reTrms

  # Response vector, design matrices
  y <- make_y_for_stan(formula, mf, family)
  x <- make_x_for_stan(formula, mf)
  z <- make_z_for_stan(formula, mf)

  # Terms
  terms <- attr(mf, "terms")
  terms <- append_predvars_attribute(terms, formula, data)

  # Binomial with >1 trials not allowed by stan_{mvmver,jm}
  is_binomial <- is.binomial(family$family)
  is_bernoulli <- is_binomial && NCOL(y$y) == 1L && all(y$y %in% 0:1)
  if (is_binomial && !is_bernoulli)
    STOP_binomial()

  # Various flags
  intercept_type <- check_intercept_type(x, family)
  has_aux <- check_for_aux(family)
  family <- append_mvmer_famlink(family, is_bernoulli)

  nlist(y, x, z, reTrms, model_frame = mf, formula, terms,
        family, intercept_type, has_aux, stub)
}
# Get the names for the Sigma var-cov matrix
#
# @param cnms The component names for the group level terms, combined
#   across all glmer submodels
# @return A character vector
get_Sigma_nms <- function(cnms) {
  nms <- names(cnms)
  Sigma_nms <- lapply(cnms, FUN = function(grp) {
    nm <- outer(grp, grp, FUN = paste, sep = ",")
    nm[lower.tri(nm, diag = TRUE)]
  })
  for (j in seq_along(Sigma_nms)) {
    Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
  }
  paste0("Sigma[", unlist(Sigma_nms), "]")
}

# Construct a list with information about the event submodel
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
handle_e_mod2 <- function(formula, data, meta, j) {

  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function")
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")

  id_var      <- meta$id_var
  id_list     <- meta$id_list
  qnodes      <- meta$qnodes[j]
  basehaz     <- meta$basehaz[[j]]
  basehaz_ops <- meta$basehaz_ops

  # parse formula, create model data & frame
  #formula   <- parse_formula(formula, data)
  formula2  <- addto_formula(formula$formula, id_var) # includes id_var
  data      <- make_model_data (formula2, data)       # row subsetting etc.
  mf_stuff  <- make_model_frame(formula2, data)       # returns Surv object

  mf <- mf_stuff$mf # model frame
  mt <- mf_stuff$mt # model terms

  mf[[id_var]] <- promote_to_factor(mf[[id_var]]) # same as lme4
  ids <- factor(mf[[id_var]])

  # error checks for the id variable
 # validate_jm_ids(y_ids = id_list, e_ids = ids)

  # entry and exit times for each row of data
  t_beg <- make_t(mf, type = "beg") # entry time
  t_end <- make_t(mf, type = "end") # exit  time
  t_upp <- make_t(mf, type = "upp") # upper time for interval censoring

  # event indicator for each row of data
  d <- make_d(mf)

  event <- as.logical(d == 1)
  rcens <- as.logical(d == 0)
  lcens <- as.logical(d == 2)
  icens <- as.logical(d == 3)

  if (any(d < 0 || d > 3))
    stop2("Invalid status indicator in Surv object.")

  # delayed entry indicator for each row of data
  delayed  <- as.logical(!t_beg == 0)

  # time variables for stan
  t_event <- t_end[event]   # exact event time
  t_lcens <- t_end[lcens]   # left  censoring time
  t_rcens <- t_end[rcens]   # right censoring time
  t_icenl <- t_end[icens]   # lower limit of interval censoring time
  t_icenu <- t_upp[icens]   # upper limit of interval censoring time
  t_delay <- t_beg[delayed] # delayed entry time

  # entry and exit times for each individual
  t_tmp <- t_end; t_tmp[icens] <- t_upp[icens]
  entrytime <- tapply(t_beg, ids, min)
  eventtime <- tapply(t_tmp, ids, max)
  status    <- tapply(d,     ids, max)

  # dimensions
  nevent <- sum(event)
  nrcens <- sum(rcens)
  nlcens <- sum(lcens)
  nicens <- sum(icens)
  ndelay <- sum(delayed)

  # baseline hazard
  ok_basehaz <- c("weibull", "bs", "piecewise")
  ok_basehaz_ops <- get_ok_basehaz_ops(basehaz)
  basehaz <- handle_basehaz(basehaz        = basehaz,
                            basehaz_ops    = basehaz_ops,
                            ok_basehaz     = ok_basehaz,
                            ok_basehaz_ops = ok_basehaz_ops,
                            times          = t_end,
                            status         = event,
                            min_t          = min(t_beg),
                            max_t          = max(c(t_end,t_upp), na.rm = TRUE))
  nvars <- basehaz$nvars # number of basehaz aux parameters

  # flag if intercept is required for baseline hazard
  has_intercept <- ai(has_intercept(basehaz))

  # standardised weights and nodes for quadrature
  qq <- get_quadpoints(nodes = qnodes)
  qp <- qq$points
  qw <- qq$weights

  # quadrature points & weights, evaluated for each row of data
  qpts_event <- uapply(qp, unstandardise_qpts, 0, t_event)
  qpts_lcens <- uapply(qp, unstandardise_qpts, 0, t_lcens)
  qpts_rcens <- uapply(qp, unstandardise_qpts, 0, t_rcens)
  qpts_icenl <- uapply(qp, unstandardise_qpts, 0, t_icenl)
  qpts_icenu <- uapply(qp, unstandardise_qpts, 0, t_icenu)
  qpts_delay <- uapply(qp, unstandardise_qpts, 0, t_delay)

  qwts_event <- uapply(qw, unstandardise_qwts, 0, t_event)
  qwts_lcens <- uapply(qw, unstandardise_qwts, 0, t_lcens)
  qwts_rcens <- uapply(qw, unstandardise_qwts, 0, t_rcens)
  qwts_icenl <- uapply(qw, unstandardise_qwts, 0, t_icenl)
  qwts_icenu <- uapply(qw, unstandardise_qwts, 0, t_icenu)
  qwts_delay <- uapply(qw, unstandardise_qwts, 0, t_delay)

  eids_event <- ids[event]
  qids_event <- rep(ids[event],   times = qnodes)
  qids_lcens <- rep(ids[lcens],   times = qnodes)
  qids_rcens <- rep(ids[rcens],   times = qnodes)
  qids_icens <- rep(ids[icens],   times = qnodes)
  qids_delay <- rep(ids[delayed], times = qnodes)

  # times at events and all quadrature points
  cids <- c(eids_event,
            qids_event,
            qids_lcens,
            qids_rcens,
            qids_icens,
            qids_icens,
            qids_delay)
  cpts_list <- list(t_event,
                    qpts_event,
                    qpts_lcens,
                    qpts_rcens,
                    qpts_icenl,
                    qpts_icenu,
                    qpts_delay)
  idx_cpts <- get_idx_array(sapply(cpts_list, length))
  cpts     <- unlist(cpts_list) # as vector for stan
  len_cpts <- length(cpts)

  # number of quadrature points
  qevent <- length(qwts_event)
  qlcens <- length(qwts_lcens)
  qrcens <- length(qwts_rcens)
  qicens <- length(qwts_icenl)
  qdelay <- length(qwts_delay)

  # basis terms for baseline hazard
  basis_cpts <- make_basis(cpts, basehaz)

  # predictor matrices
  x <- make_x(formula$tf_form, mf)$x
  x_event <- keep_rows(x, d == 1)
  x_lcens <- keep_rows(x, d == 2)
  x_rcens <- keep_rows(x, d == 0)
  x_icens <- keep_rows(x, d == 3)
  x_delay <- keep_rows(x, delayed)
  K <- ncol(x)
  x_cpts <- rbind(x_event,
                  rep_rows(x_event, times = qnodes),
                  rep_rows(x_lcens, times = qnodes),
                  rep_rows(x_rcens, times = qnodes),
                  rep_rows(x_icens, times = qnodes),
                  rep_rows(x_delay, times = qnodes))

  # fit a cox model
  if (formula$surv_type %in% c("right", "counting")) {
    mod <- survival::coxph(formula$formula, data = data, x = TRUE)
  } else if (formula$surv_type %in% c("interval", "interval2")) {
    mod <- survival::survreg(formula$formula, data = data, x = TRUE)
  } else {
    stop("Bug found: Invalid Surv type.")
  }

  # calculate mean log incidence, used as a shift in log baseline hazard
  norm_const <- log(nevent / sum(eventtime - entrytime))

  nlist(mod,
        surv_type = formula$surv_type,
        qnodes,
        basehaz,
        has_intercept,
        has_icens = as.logical(nicens),
        model_frame = mf,
        entrytime,
        eventtime,
        d, status,
        norm_const,
        t_beg,
        t_end,
        t_upp,
        t_event,
        t_rcens,
        t_lcens,
        t_icenl,
        t_icenu,
        t_delay,
        nevent,
        nlcens,
        nrcens,
        nicens,
        ndelay,
        qevent,
        qlcens,
        qrcens,
        qicens,
        qdelay,
        cids,
        cpts,
        len_cpts,
        idx_cpts,
        qwts_event,
        qwts_lcens,
        qwts_rcens,
        qwts_icenl,
        qwts_icenu,
        qwts_delay,
        eids_event,
        qids_event,
        qids_lcens,
        qids_rcens,
        qids_icens,
        qids_delay,
        x,
        x_cpts,
        basis_cpts,
        x_bar = colMeans(x),
        K)
}


# Check that the observation times for the longitudinal submodel are all
# positive and not observed after the individual's event time
#
# @param data A data frame (data for one longitudinal submodel)
# @param eventtimes A named numeric vector with the event time for each
#   individual. The vector names should be the individual ids.
# @param id_var,time_var The ID and time variable in the longitudinal data.
# @return Nothing.
validate_observation_times <-function(data, eventtimes, id_var, time_var) {
  if (!time_var %in% colnames(data))
    STOP_no_var(time_var)
  if (!id_var %in% colnames(data))
    STOP_no_var(id_var)
  if (any(data[[time_var]] < 0))
    stop2("Values for the time variable (", time_var, ") should not be negative.")
  mt  <- tapply(data[[time_var]], factor(data[[id_var]]), max) # max observation time
  nms <- names(eventtimes)                                     # patient IDs
  if (is.null(nms))
    stop2("Bug found: cannot find names in the vector of event times.")
  sel <- which(sapply(nms, FUN = function(i) mt[i] > eventtimes[i]))
  if (length(sel))
    stop2("The following individuals have observation times in the longitudinal ",
          "data that are later than their event time: ", comma(nms[sel]))
}

# Filter associated times between the longitudinal model and hazard sumbmodel
#
# @param x Long dataset
# @param t_var time variable from long dataset.
# @param t time event vector
assoc_time <- function(x, t_var, t, id_state){
  tmpdf <- data.frame(id = id_state,
                      eventtime = t)
  tmpx <- dplyr::left_join(x, tmpdf, by = "id")

  out <- tmpx[[t_var]] <= tmpx[["eventtime"]]
  out[is.na(out)] <- FALSE
  out
}

extract_id <- function(x, id_var){
  x[[id_var]]
}


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
handle_ms_mod <- function(formula, data, meta) {

  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function")
  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")

  id_var      <- meta$id_var
  id_list     <- meta$id_list
  time_start  <- meta$time_start
  qnodes      <- meta$qnodes
  basehaz     <- meta$basehaz
  basehaz_ops <- meta$basehaz_ops

  # parse formula, create model data & frame
  formula2  <- lapply(formula, function(f) addto_formula(f$formula, id_var)) # includes id_var

  mf_stuff <-  lapply(seq_along(formula2), function(i) {
    make_model_frame(formula2[[i]], data[[i]])
  })

  mf <- lapply(mf_stuff, function(m) m$mf)   # model frame
  mt <- lapply(mf_stuff, function(m) m$mt)  # model terms

  for(i in seq_along(mf)){
    mf[[i]][[id_var]] <- promote_to_factor(mf[[i]][[id_var]]) # same as lme4
  }

  ids <- lapply(mf, function(mf) factor(mf[[id_var]]) )

  # error checks for the id variable
  validate_jm_ids(y_ids = id_list, e_ids = uu(ids) )

  # entry and exit times for each row of data
  t_beg <- lapply(mf, function(m) make_t(m, type = "beg") ) # entry time
  t_end <- lapply(mf, function(m) make_t(m, type = "end") ) # exit time
  t_upp <- lapply(mf, function(m) make_t(m, type = "upp") ) # upper time for

  # ensure no event or censoring times are zero (leads to degenerate
  # estimate for log hazard for most baseline hazards, due to log(0))

  for(i in seq_along(t_end ) ){
    check1 <- any(t_end[[i]] <= 0, na.rm = TRUE)
    if (check1)
      stop2("All event and censoring times must be greater than 0.")
  }

  # event indicator for each row of data
  d <- lapply(mf, make_d)

  for(i in seq_along(status)){
    if (any(status[[i]] < 0 || status[[i]] > 3))
      stop2("Invalid status indicator in formula.")
  }

  event <- lapply(d, function(d) as.logical(d == 1) )
  rcens <- lapply(d, function(d) as.logical(d == 0) )
  lcens <- lapply(d, function(d) as.logical(d == 2) )
  icens <- lapply(d, function(d) as.logical(d == 3) )

  # delayed entry indicator for each row of data
  delayed  <- lapply(t_beg, function(t_beg) as.logical(!t_beg == 0) )

  # time variables for stan
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


  # entry and exit times for each individual
  t_tmp <- t_end
  entrytime <- list(); eventtime <- list();
  for(i in seq_along(t_tmp)){
    t_tmp[[i]][icens[[i]]] <- t_upp[[i]][icens[[i]]]
    entrytime[[i]] <- tapply(t_beg[[i]], ids[[i]], min)
    eventtime[[i]] <- tapply(t_tmp[[i]], ids[[i]], max)
    status[[i]]    <- tapply(d[[i]],     ids[[i]], max)
  }

  # dimensions
  nevent <- lapply(status, function(s)  sum(s == 1))
  nrcens <- lapply(status,  function(s) sum(s == 0))
  nlcens <- lapply(status, function(s) sum(s == 2))
  nicens <- lapply(status, function(s) sum(s == 3))
  ndelay <- lapply(delayed, function(d) sum(d))

  # baseline hazard
  ok_basehaz <- c("weibull", "bs", "piecewise")
  ok_basehaz_ops <- lapply(basehaz, function(b) get_ok_basehaz_ops(b))

  basehaz <- lapply(seq_along(basehaz), function(b)
    SW( handle_basehaz_surv(basehaz = basehaz[[b]],
                            basehaz_ops = basehaz_ops[[b]],
                            ok_basehaz = ok_basehaz,
                            ok_basehaz_ops = ok_basehaz_ops[[b]],
                            times = t_end[[b]],
                            status = event[[b]],
                            min_t = min(t_beg[[b]]),
                            max_t = max(c(t_end[[b]], t_upp[[b]]), na.rm = TRUE ) )
    ))

  nvars <- lapply(basehaz, function(b) b$nvars)  # number of basehaz aux parameters
  # flag if intercept is required for baseline hazard
  has_intercept   <- lapply(basehaz, function(b) ai(has_intercept(b)) )

  nvars <- basehaz$nvars # number of basehaz aux parameters

  # flag if intercept is required for baseline hazard
  has_intercept <- ai(has_intercept(basehaz))

  # standardised weights and nodes for quadrature
  qq <- lapply(qnodes, get_quadpoints)
  qp <- lapply(qq, function(qq) qq$points)
  qw <- lapply(qq, function(qq) qq$weights)

  # quadrature points & weights, evaluated for each row of data
  qpts_event <- list()
  qpts_lcens <- list()
  qpts_rcens <- list()
  qpts_icenl <- list()
  qpts_icenu <- list()
  qpts_delay <- list()

  qwts_event <- list()
  qwts_lcens <- list()
  qwts_rcens <- list()
  qwts_icenl <- list()
  qwts_icenu <- list()
  qwts_delay <- list()

  for(i in seq_along(qq)){
    qp_i = qp[[i]]
    qw_i = qw[[i]]
    qpts_event[[i]] <- uapply(qp_i,unstandardise_qpts, 0, t_event[[i]])
    qpts_lcens[[i]] <- uapply(qp_i,unstandardise_qpts, 0, t_lcens[[i]])
    qpts_rcens[[i]] <- uapply(qp_i,unstandardise_qpts, 0, t_rcens[[i]])
    qpts_icenl[[i]] <- uapply(qp_i,unstandardise_qpts, 0, t_icenl[[i]])
    qpts_icenu[[i]] <- uapply(qp_i,unstandardise_qpts, 0, t_icenu[[i]])
    qpts_delay[[i]] <- uapply(qp_i,unstandardise_qpts, 0, t_delay[[i]])

    qwts_event[[i]] <- uapply(qw_i, unstandardise_qwts, 0, t_event[[i]])
    qwts_lcens[[i]] <- uapply(qw_i, unstandardise_qwts, 0, t_lcens[[i]])
    qwts_rcens[[i]] <- uapply(qw_i, unstandardise_qwts, 0, t_rcens[[i]])
    qwts_icenl[[i]] <- uapply(qw_i, unstandardise_qwts, 0, t_icenl[[i]])
    qwts_icenu[[i]] <- uapply(qw_i, unstandardise_qwts, 0, t_icenu[[i]])
    qwts_delay[[i]] <- uapply(qw_i, unstandardise_qwts, 0, t_delay[[i]])
  }


  eids_event <- ids[event]
  qids_event <- rep(ids[event],   times = qnodes)
  qids_lcens <- rep(ids[lcens],   times = qnodes)
  qids_rcens <- rep(ids[rcens],   times = qnodes)
  qids_icens <- rep(ids[icens],   times = qnodes)
  qids_delay <- rep(ids[delayed], times = qnodes)

  # times at events and all quadrature points
  cids <- c(eids_event,
            qids_event,
            qids_lcens,
            qids_rcens,
            qids_icens,
            qids_icens,
            qids_delay)
  cpts_list <- list(t_event,
                    qpts_event,
                    qpts_lcens,
                    qpts_rcens,
                    qpts_icenl,
                    qpts_icenu,
                    qpts_delay)
  idx_cpts <- get_idx_array(sapply(cpts_list, length))
  cpts     <- unlist(cpts_list) # as vector for stan
  len_cpts <- length(cpts)

  # number of quadrature points
  qevent <- length(qwts_event)
  qlcens <- length(qwts_lcens)
  qrcens <- length(qwts_rcens)
  qicens <- length(qwts_icenl)
  qdelay <- length(qwts_delay)

  # basis terms for baseline hazard
  basis_cpts <- make_basis(cpts, basehaz)

  # predictor matrices
  x <- make_x(formula$tf_form, mf)$x
  x_event <- keep_rows(x, d == 1)
  x_lcens <- keep_rows(x, d == 2)
  x_rcens <- keep_rows(x, d == 0)
  x_icens <- keep_rows(x, d == 3)
  x_delay <- keep_rows(x, delayed)
  K <- ncol(x)
  x_cpts <- rbind(x_event,
                  rep_rows(x_event, times = qnodes),
                  rep_rows(x_lcens, times = qnodes),
                  rep_rows(x_rcens, times = qnodes),
                  rep_rows(x_icens, times = qnodes),
                  rep_rows(x_delay, times = qnodes))

  # fit a cox model
  if (formula$surv_type %in% c("right", "counting")) {
    mod <- survival::coxph(formula$formula, data = data, x = TRUE)
  } else if (formula$surv_type %in% c("interval", "interval2")) {
    mod <- survival::survreg(formula$formula, data = data, x = TRUE)
  } else {
    stop("Bug found: Invalid Surv type.")
  }

  # calculate mean log incidence, used as a shift in log baseline hazard
  norm_const <- log(nevent / sum(eventtime - entrytime))

  nlist(mod,
        surv_type = formula$surv_type,
        qnodes,
        basehaz,
        has_intercept,
        has_icens = as.logical(nicens),
        model_frame = mf,
        entrytime,
        eventtime,
        d, status,
        norm_const,
        t_beg,
        t_end,
        t_upp,
        t_event,
        t_rcens,
        t_lcens,
        t_icenl,
        t_icenu,
        t_delay,
        time_start,
        nevent,
        nlcens,
        nrcens,
        nicens,
        ndelay,
        qevent,
        qlcens,
        qrcens,
        qicens,
        qdelay,
        cids,
        cpts,
        len_cpts,
        idx_cpts,
        qwts_event,
        qwts_lcens,
        qwts_rcens,
        qwts_icenl,
        qwts_icenu,
        qwts_delay,
        eids_event,
        qids_event,
        qids_lcens,
        qids_rcens,
        qids_icens,
        qids_delay,
        x,
        x_cpts,
        basis_cpts,
        x_bar = colMeans(x),
        K)
}

# Construct a list with information about the baseline hazard
#
# @param basehaz A string specifying the type of baseline hazard
# @param basehaz_ops A named list with elements df, knots
# @param ok_basehaz A list of admissible baseline hazards
# @param times A numeric vector with eventtimes for each individual
# @param status A numeric vector with event indicators for each individual
# @param min_t Scalar, the minimum entry time across all individuals
# @param max_t Scalar, the maximum event or censoring time across all individuals
# @return A named list with the following elements:
#   type: integer specifying the type of baseline hazard, 1L = weibull,
#     2L = b-splines, 3L = piecewise.
#   type_name: character string specifying the type of baseline hazard.
#   user_df: integer specifying the input to the df argument
#   df: integer specifying the number of parameters to use for the
#     baseline hazard.
#   knots: the knot locations for the baseline hazard.
#   bs_basis: The basis terms for the B-splines. This is passed to Stan
#     as the "model matrix" for the baseline hazard. It is also used in
#     post-estimation when evaluating the baseline hazard for posterior
#     predictions since it contains information about the knot locations
#     for the baseline hazard (this is implemented via splines::predict.bs).
handle_basehaz <- function(basehaz,
                           basehaz_ops,
                           ok_basehaz     = c("weibull", "bs", "piecewise"),
                           ok_basehaz_ops = c("df", "knots"),
                           times,
                           status,
                           min_t, max_t) {

  if (!basehaz %in% ok_basehaz)
    stop2("'basehaz' should be one of: ", comma(ok_basehaz))

  if (!all(names(basehaz_ops) %in% ok_basehaz_ops))
    stop2("'basehaz_ops' can only include: ", comma(ok_basehaz_ops))

  if (basehaz == "exp") {

    bknots <- NULL # boundary knot locations
    iknots <- NULL # internal knot locations
    basis  <- NULL # spline basis
    nvars  <- 0L   # number of aux parameters, none

  } else if (basehaz == "gompertz") {

    bknots <- NULL # boundary knot locations
    iknots <- NULL # internal knot locations
    basis  <- NULL # spline basis
    nvars  <- 1L   # number of aux parameters, Gompertz scale

  } else if (basehaz == "weibull") {

    bknots <- NULL # boundary knot locations
    iknots <- NULL # internal knot locations
    basis  <- NULL # spline basis
    nvars  <- 1L   # number of aux parameters, Weibull shape

  } else if (basehaz == "bs") {

    df    <- basehaz_ops$df
    knots <- basehaz_ops$knots

    if (!is.null(df) && !is.null(knots))
      stop2("Cannot specify both 'df' and 'knots' for the baseline hazard.")

    if (is.null(df))
      df <- 5L # default df for B-splines, assuming no intercept
    # NB this is ignored if the user specified knots

    tt <- times[status == 1] # uncensored event times
    if (is.null(knots) && !length(tt)) {
      warning2("No observed events found in the data. Censoring times will ",
               "be used to evaluate default knot locations for splines.")
      tt <- times
    }

    bknots <- c(min_t, max_t)
    iknots <- get_iknots(tt, df = df, iknots = knots)
    basis  <- get_basis(tt, iknots = iknots, bknots = bknots, type = "bs")
    nvars  <- ncol(basis)  # number of aux parameters, basis terms

  } else if (basehaz == "ms") {

    df    <- basehaz_ops$df
    knots <- basehaz_ops$knots

    if (!is.null(df) && !is.null(knots)) {
      stop2("Cannot specify both 'df' and 'knots' for the baseline hazard.")
    }

    tt <- times[status == 1] # uncensored event times
    if (is.null(df)) {
      df <- 5L # default df for B-splines, assuming no intercept
      # NB this is ignored if the user specified knots
    }

    tt <- times[status == 1] # uncensored event times
    if (is.null(knots) && !length(tt)) {
      warning2("No observed events found in the data. Censoring times will ",
               "be used to evaluate default knot locations for splines.")
      tt <- times
    }

    bknots <- c(min_t, max_t)
    iknots <- get_iknots(tt, df = df, iknots = knots)
    basis  <- get_basis(tt, iknots = iknots, bknots = bknots, type = "ms")
    nvars  <- ncol(basis)  # number of aux parameters, basis terms

  } else if (basehaz == "piecewise") {

    df    <- basehaz_ops$df
    knots <- basehaz_ops$knots

    if (!is.null(df) && !is.null(knots)) {
      stop2("Cannot specify both 'df' and 'knots' for the baseline hazard.")
    }

    if (is.null(df)) {
      df <- 6L # default number of segments for piecewise constant
      # NB this is ignored if the user specified knots
    }

    if (is.null(knots) && !length(tt)) {
      warning2("No observed events found in the data. Censoring times will ",
               "be used to evaluate default knot locations for piecewise basehaz.")
      tt <- times
    }

    bknots <- c(min_t, max_t)
    iknots <- get_iknots(tt, df = df, iknots = knots)
    basis  <- NULL               # spline basis
    nvars  <- length(iknots) + 1 # number of aux parameters, dummy indicators

  }

  nlist(type_name = basehaz,
        type = basehaz_for_stan(basehaz),
        nvars,
        iknots,
        bknots,
        basis,
        df = nvars,
        user_df = nvars,
        knots = if (basehaz == "bs") iknots else c(bknots[1], iknots, bknots[2]),
        bs_basis = basis)
}

# Append a family object with numeric family and link information used by Stan
#
# @param family The existing family object
# @param is_bernoulli Logical specifying whether the family should be bernoulli
# @return A family object with two appended elements:
#   mvmer_family: an integer telling Stan which family
#   mvmer_link: an integer telling Stan which link function (varies by family!)
append_mvmer_famlink <- function(family, is_bernoulli = FALSE) {
  famname <- family$family
  family$mvmer_family <- switch(
    famname,
    gaussian = 1L,
    Gamma = 2L,
    inverse.gaussian = 3L,
    binomial = 5L, # bernoulli = 4L changed later
    poisson = 6L,
    "neg_binomial_2" = 7L)
  if (is_bernoulli)
    family$mvmer_family <- 4L
  supported_links <- supported_glm_links(famname)
  link <- which(supported_links == family$link)
  family$mvmer_link <- link
  return(family)
}


# Function to check if the submodel should include a auxiliary term
#
# @param family A family object
# @return A logical specify whether the submodel includes a auxiliary term
check_for_aux <- function(family) {
  !(family$family %in% c("binomial", "poisson"))
}


#--------------- Functions related to longitudinal submodel

# Return the response vector for passing to Stan
#
# @param formula The model formula
# @param model_frame The model frame
# @param family A family object
# @return A named list with the following elements:
#   y: the response vector
#   real: the response vector if real, else numeric(0)
#   integer: the response vector if integer, else integer(0)
#   resp_type: 1L if response is real, 2L is response is integer
make_y_for_stan <- function(formula, model_frame, family) {
  y <- as.vector(model.response(model_frame))
  y <- validate_glm_outcome_support(y, family)
  resp_type <- if (check_response_real(family)) 1L else 2L
  real    <- if (resp_type == 1L) y else numeric(0)
  integer <- if (resp_type == 2L) y else integer(0)
  nlist(y, real, integer, resp_type)
}

# Return the design matrix for passing to Stan
#
# @param formula The model formula.
# @param model_frame The model frame.
# @return A named list with the following elements:
#   x: the fe model matrix, not centred and may have intercept.
#   xtemp: fe model matrix, centred and no intercept.
#   x_form: the formula for the fe model matrix.
#   x_bar: the column means of the model matrix.
#   has_intercept: logical for whether the submodel has an intercept
#   N,K: number of rows (observations) and columns (predictors) in the
#     fixed effects model matrix
make_x_for_stan <- function(formula, model_frame) {
  x_form <- lme4::nobars(formula)
  x <- model.matrix(x_form, model_frame)
  has_intercept <- check_for_intercept(x, logical = TRUE)
  xtemp <- drop_intercept(x)
  x_bar <- colMeans(xtemp)
  xtemp <- sweep(xtemp, 2, x_bar, FUN = "-")
  # identify any column of x with < 2 unique values (empty interaction levels)
  sel <- (2 > apply(xtemp, 2L, function(x) length(unique(x))))
  if (any(sel))
    stop2("Cannot deal with empty interaction levels found in columns: ",
          paste(colnames(xtemp)[sel], collapse = ", "))
  nlist(x, xtemp, x_form, x_bar, has_intercept, N = NROW(xtemp), K = NCOL(xtemp))
}

# Split the random effects part of a model formula into
#   - the formula part (ie. the formula on the LHS of "|"), and
#   - the name of the grouping factor (ie. the variable on the RHS of "|")
#
# @param x Random effects part of a model formula, as returned by lme4::findbars
# @return A named list with the following elements:
#   re_form: a formula specifying the random effects structure
#   group_var: the name of the grouping factor
split_at_bars <- function(x) {
  terms <- strsplit(deparse(x, 500), "\\s\\|\\s")[[1L]]
  if (!length(terms) == 2L)
    stop2("Could not parse the random effects formula.")
  re_form <- formula(paste("~", terms[[1L]]))
  group_var <- terms[[2L]]
  nlist(re_form, group_var)
}

# Return design matrices for the group level terms for passing to Stan
#
# @param formula The model formula
# @param model_frame The model frame
# @return A named list with the following elements:
#   z: a list with each element containing the random effects model
#     matrix for one grouping factor.
#   z_forms: a list with each element containing the model formula for
#     one grouping factor.
#   group_vars: a character vector with the name of each of the
#     grouping factors
#   group_cnms: a list with each element containing the names of the
#     group level parameters for one grouping factor
#   group_list: a list with each element containing the vector of group
#     IDs for the rows of z
#   nvars: a vector with the number of group level parameters for each
#     grouping factor
#   ngrps: a vector with the number of groups for each grouping factor
make_z_for_stan <- function(formula, model_frame) {
  bars <- lme4::findbars(formula)
  if (length(bars) > 2L)
    stop2("A maximum of 2 grouping factors are allowed.")
  z_parts <- lapply(bars, split_at_bars)
  z_forms <- fetch(z_parts, "re_form")
  z <- lapply(z_forms, model.matrix, model_frame)
  group_cnms <- lapply(z, colnames)
  group_vars <- fetch(z_parts, "group_var")
  group_list <- lapply(group_vars, function(x) factor(model_frame[[x]]))
  nvars <- lapply(group_cnms, length)
  ngrps <- lapply(group_list, n_distinct)
  names(z) <- names(z_forms) <- names(group_cnms) <-
    names(group_list) <- names(nvars) <- names(ngrps) <- group_vars
  nlist(z, z_forms, group_vars, group_cnms, group_list, nvars, ngrps)
}


# Take the model frame terms object and append with attributes
# that provide the predvars for the fixed and random effects
# parts, based on the model formula and data
#
# @param terms The existing model frame terms object
# @param formula The formula that was used to build the model frame
#   (but prior to having called lme4::subbars on it!)
# @param data The data frame that was used to build the model frame
# @return A terms object with predvars.fixed and predvars.random as
#   additional attributes
append_predvars_attribute <- function(terms, formula, data) {
  fe_form <- lme4::nobars(formula)
  re_form <- lme4::subbars(justRE(formula, response = TRUE))
  fe_frame <- stats::model.frame(fe_form, data)
  re_frame <- stats::model.frame(re_form, data)
  fe_terms <- attr(fe_frame, "terms")
  re_terms <- attr(re_frame, "terms")
  fe_predvars <- attr(fe_terms, "predvars")
  re_predvars <- attr(re_terms, "predvars")
  attr(terms, "predvars.fixed")  <- attr(fe_terms, "predvars")
  attr(terms, "predvars.random") <- attr(re_terms, "predvars")
  terms
}


# Return info on the required type of intercept
#
# @param X The model matrix
# @param family A family object
# @return A named list with the following elements:
#   type: character string specifying the type of bounds to use
#     for the intercept.
#   number: an integer specifying the type of bounds to use
#     for the intercept where 0L = no intercept, 1L = no bounds
#     on intercept, 2L = lower bound, 3L = upper bound.
check_intercept_type <- function(X, family) {
  fam <- family$family
  link <- family$link
  if (!X$has_intercept) { # no intercept
    type <- "none"
    needs_intercept <-
      (!is.gaussian(fam) && link == "identity") ||
      (is.gamma(fam) && link == "inverse") ||
      (is.binomial(fam) && link == "log")
    if (needs_intercept)
      stop2("To use the specified combination of family and link (", fam,
            ", ", link, ") the model must have an intercept.")
  } else if (fam == "binomial" && link == "log") { # binomial, log
    type <- "upper_bound"
  } else if (fam == "binomial") { # binomial, !log
    type <- "no_bound"
  } else if (link == "log") { # gamma/inv-gaus/poisson/nb, log
    type <- "no_bound"
  } else if (fam == "gaussian") { # gaussian, !log
    type <- "no_bound"
  } else { # gamma/inv-gaus/poisson/nb, !log
    type <- "lower_bound"
  }
  number <- switch(type, none = 0L, no_bound = 1L,
                   lower_bound = 2L, upper_bound = 3L)
  nlist(type, number)
}

justRE <- function(f, response = FALSE) {
  response <- if (response && length(f) == 3) f[[2]] else NULL
  reformulate(paste0("(", vapply(lme4::findbars(f),
                                 function(x) paste(deparse(x, 500L),
                                                   collapse = " "),
                                 ""), ")"),
              response = response)
}


# Function to return a single cnms object for all longitudinal submodels
#
# @param x A list, with each element being a cnms object returned by (g)lmer
get_common_cnms <- function(y_mod, stub = "Long") {
  y_cnms <- fetch(y_mod, "z", "group_cnms")
  nms <- lapply(y_cnms, names)
  unique_nms <- unique(unlist(nms))
  cnms <- lapply(seq_along(unique_nms), function(i) {
    nm <- unique_nms[i]
    unlist(lapply(1:length(y_cnms), function(m)
      if (nm %in% nms[[m]]) paste0(stub, m, "|", y_cnms[[m]][[nm]])))
  })
  names(cnms) <- unique_nms
  if (length(cnms) > 2L)
    stop("A maximum of 2 grouping factors are allowed.")
  cnms
}

# Function to return a single list with the factor levels for each
# grouping factor, but collapsed across all longitudinal submodels
#
# @param x A list containing the flist object for each of the submodels
get_common_flevels <- function(y_mod) {
  y_flevels <- fetch(y_mod, "z", "group_list")
  nms <- lapply(y_flevels, names)
  unique_nms <- unique(unlist(nms))
  flevels <- lapply(seq_along(unique_nms), function(i) {
    nm <- unique_nms[i]
    flevels_nm <- lapply(1:length(y_flevels), function(m)
      if (nm %in% nms[[m]]) levels(y_flevels[[m]][[nm]]))
    flevels_nm <- rm_null(unique(flevels_nm))
    if (length(flevels_nm) > 1L)
      stop2("The group factor levels must be the same for all submodels.")
    flevels_nm[[1L]]
  })
  names(flevels) <- unique_nms
  flevels
}

# Check the id_var argument is valid and is included appropriately in the
# formulas for each of the longitudinal submodels
#
# @param id_var The character string that the user specified for the id_var
#   argument -- will have been set to NULL if the argument was missing.
# @param y_cnms A list of length M with the cnms for each longitudinal submodel
# @param y_flist A list of length M with the flist for each longitudinal submodel
# @return Returns the character string corresponding to the appropriate id_var.
#   This will either be the user specified id_var argument or the only grouping
#   factor.
check_id_var <- function(y_mod, id_var) {

  y_cnms <- fetch(y_mod, "z", "group_cnms")

  len_cnms <- sapply(y_cnms, length) # num grouping factors in each submodel

  if (any(len_cnms > 1L)) { # multiple grouping factors

    if (is.null(id_var))
      stop2("With more than one grouping factor 'id_var' must be specified.")
    lapply(y_cnms, function(x)  if (!(id_var %in% names(x)))
      stop2("'id_var' must be a grouping factor in each longitudinal submodel."))
    return(id_var)

  } else { # only one grouping factor (assumed to be subject ID)

    only_cnm <- unique(sapply(y_cnms, names))
    if (length(only_cnm) > 1L)
      stop2("The grouping factor (ie, subject ID variable) is not the ",
            "same in all longitudinal submodels.")
    if (not.null(id_var) && !identical(id_var, only_cnm))
      warning2("The user specified 'id_var' (", paste(id_var),
               ") and the assumed ID variable based on the single ",
               "grouping factor (", paste(only_cnm), ") are not the same; ",
               "'id_var' will be ignored")
    return(only_cnm)
  }
}

# Check the factor list corresponding to subject ID is the same in each
# of the longitudinal submodels
#
# @param id_var The name of the ID variable
# @param y_flist A list containing the flist objects returned for each
#   separate longitudinal submodel
# @return A vector of factor levels corresponding to the IDs appearing
#   in the longitudinal submodels
check_id_list <- function(y_mod, id_var) {
  y_flist <- fetch(y_mod, "z", "group_list")
  id_list <- unique(lapply(y_flist, function(x) levels(x[[id_var]])))
  if (length(id_list) > 1L)
    stop2("The subject IDs are not the same in all longitudinal submodels.")
  unlist(id_list)
}

# Check the weights argument for stan_jm
#
# @param weights The data frame passed via the weights argument
# @param id_var The name of the ID variable
check_weights <- function(weights, id_var) {

  if (is.null(weights))
    return(weights)

  # Check weights are an appropriate data frame
  if ((!is.data.frame(weights)) || (!ncol(weights) == 2))
    stop("'weights' argument should be a data frame with two columns: the first ",
         "containing patient IDs, the second containing their corresponding ",
         "weights.", call. = FALSE)
  if (!id_var %in% colnames(weights))
    stop("The data frame supplied in the 'weights' argument should have a ",
         "column named ", id_var, call. = FALSE)
  weight_var <- setdiff(colnames(weights), id_var)

  # Check weights are positive and numeric
  wts <- weights[[weight_var]]
  if (!is.numeric(wts))
    stop("The weights supplied must be numeric.", call. = FALSE)
  if (any(wts < 0))
    stop("Negative weights are not allowed.", call. = FALSE)

  # Check only one weight per ID
  n_weights_per_id <- tapply(weights[[weight_var]], weights[[id_var]], length)
  if (!all(n_weights_per_id == 1L))
    stop("The data frame supplied in the 'weights' argument should only have ",
         "one row (ie, one weight) per patient ID.", call. = FALSE)
}

# Return the vector of prior weights for one of the submodels
#
# @param mod_stuff A named list with elements: y, flist, ord
# @param weights The data frame passed via the weights argument
# @param id_var The name of the ID variable
handle_weights <- function(mod_stuff, weights, id_var) {

  is_glmod <- (is.null(mod_stuff$eventtime))

  # No weights provided by user
  if (is.null(weights)) {
    len <- if (is_glmod) length(mod_stuff$Y$Y) else 0
    return(rep(0.0, len))
  }

  # Check for IDs with no weight supplied
  weights[[id_var]] <- factor(weights[[id_var]])
  ids <- if (is_glmod) mod_stuff$Z$group_list[[id_var]] else factor(mod_stuff$id_list)
  sel <- which(!ids %in% weights[[id_var]])
  if (length(sel)) {
    if (length(sel) > 30L) sel <- sel[1:30]
    stop(paste0("The following patient IDs are used in fitting the model, but ",
                "do not have weights supplied via the 'weights' argument: ",
                paste(ids[sel], collapse = ", ")), call. = FALSE)
  }

  # Obtain length and ordering of weights vector using flist
  wts_df  <- merge(data.frame(id = ids), weights, by.x = "id", by.y = id_var, sort = FALSE)
  wts_var <- setdiff(colnames(weights), id_var)
  wts     <- wts_df[[wts_var]]

  wts
}


# Deal with covariance prior
#
# @param prior A list
# @param cnms A list of lists, with names of the group specific
#   terms for each grouping factor
# @param ok_dists A list of admissible distributions
handle_cov_prior <- function(prior, cnms, ok_dists = nlist("decov", "lkj")) {
  if (!is.list(prior))
    stop(sQuote(deparse(substitute(prior))), " should be a named list")
  t <- length(unique(cnms)) # num grouping factors
  p <- sapply(cnms, length) # num terms for each grouping factor
  prior_dist_name <- prior$dist
  if (!prior_dist_name %in% unlist(ok_dists)) {
    stop("The prior distribution should be one of ",
         paste(names(ok_dists), collapse = ", "))
  } else if (prior_dist_name == "decov") {
    prior_shape <- as.array(maybe_broadcast(prior$shape, t))
    prior_scale <- as.array(maybe_broadcast(prior$scale, t))
    prior_concentration <-
      as.array(maybe_broadcast(prior$concentration, sum(p[p > 1])))
    prior_regularization <-
      as.array(maybe_broadcast(prior$regularization, sum(p > 1)))
    prior_df <- NULL
  } else if (prior_dist_name == "lkj") {
    prior_shape <- NULL
    prior_scale <- as.array(maybe_broadcast(prior$scale, sum(p)))
    prior_concentration <- NULL
    prior_regularization <-
      as.array(maybe_broadcast(prior$regularization, sum(p > 1)))
    prior_df <- as.array(maybe_broadcast(prior$df, sum(p)))
  }
  prior_dist <- switch(prior_dist_name, decov = 1L, lkj = 2L)

  nlist(prior_dist_name, prior_dist, prior_shape, prior_scale,
        prior_concentration, prior_regularization, prior_df, t, p,
        prior_autoscale = isTRUE(prior$autoscale))
}

# Seperate the information about the covariance prior into a list
# of lists. At the top level of the returned list the elements
# correpond to each of the grouping factors, and on the second level
# of the returned list the elements correpsond to the separate glmer
# submodels. This separation is required for autoscaling the priors
# on the sds of group level effects, since these are autoscaled based
# on the separate Z matrices (design matrices for the random effects).
#
# @param prior_stuff The named list returned by handle_cov_prior
# @param cnms The component names for group level terms, combined across
#   all glmer submodels
# @param submodel_cnms The component names for the group level terms,
#   separately for each glmer submodel (stored as a list of length M)
# @return A list with each element containing the covariance prior
#   information for one grouping factor
split_cov_prior <- function(prior_stuff, cnms, submodel_cnms) {
  if (!prior_stuff$prior_dist_name == "lkj") {
    return(prior_stuff) # nothing to be done for decov prior
  } else {
    M <- length(submodel_cnms) # number of submodels
    cnms_nms <- names(cnms) # names of grouping factors
    mark <- 0
    new_prior_stuff <- list()
    for (nm in cnms_nms) {
      for (m in 1:M) {
        len <- length(submodel_cnms[[m]][[nm]])
        new_prior_stuff[[nm]][[m]] <- prior_stuff
        if (len) {
          # submodel 'm' has group level terms for group factor 'nm'
          beg <- mark + 1; end <- mark + len
          new_prior_stuff[[nm]][[m]]$prior_scale <- prior_stuff$prior_scale[beg:end]
          new_prior_stuff[[nm]][[m]]$prior_df <- prior_stuff$prior_df[beg:end]
          mark <- mark + len
        } else {
          new_prior_stuff[[nm]][[m]]$prior_scale <- NULL
          new_prior_stuff[[nm]][[m]]$prior_df <- NULL
          new_prior_stuff[[nm]][[m]]$prior_regularization <- NULL
        }
      }
    }
  }
  new_prior_stuff
}

# From a vector of length M giving the number of elements (for example number
# of parameters or observations) for each submodel, create an indexing array
# of dimension M * 2, where column 1 is the beginning index and 2 is the end index
#
# @param x A numeric vector
# @return A length(x) * 2 array
get_idx_array <- function(x) {
  aa(do.call("rbind", lapply(1:length(x), function(i) {
    idx_beg <- ifelse(x[i] > 0L, sum(x[0:(i-1)]) + 1, 0L)
    idx_end <- ifelse(x[i] > 0L, sum(x[0:i]),         0L)
    c(idx_beg, idx_end)
  })))
}

# Check that the ids in the longitudinal and survival models match
validate_jm_ids <- function(y_ids, e_ids) {
  if (!identical(y_ids, levels(factor(e_ids))))
    stop2("The patient IDs (levels of the grouping factor) included ",
          "in the longitudinal and event submodels do not match")
  if (is.unsorted(factor(e_ids)))
    stop2("'dataEvent' needs to be sorted by the subject ID/grouping variable.")
  if (!identical(length(y_ids), length(e_ids)))
    stop2("The number of patients differs between the longitudinal and ",
          "event submodels. Perhaps you intended to use 'start/stop' notation ",
          "for the Surv() object.")
}


# Function to return standardised GK quadrature points and weights
#
# @param nodes The required number of quadrature nodes
# @return A list with two named elements (points and weights) each
#   of which is a numeric vector with length equal to the number of
#   quadrature nodes
get_quadpoints <- function(nodes = 15) {
  if (!is.numeric(nodes) || (length(nodes) > 1L)) {
    stop("'qnodes' should be a numeric vector of length 1.")
  } else if (nodes == 15) {
    list(
      points = c(
        -0.991455371120812639207,
        -0.949107912342758524526,
        -0.86486442335976907279,
        -0.7415311855993944398639,
        -0.5860872354676911302941,
        -0.4058451513773971669066,
        -0.2077849550078984676007,
        0,
        0.2077849550078984676007,
        0.405845151377397166907,
        0.5860872354676911302941,
        0.741531185599394439864,
        0.86486442335976907279,
        0.9491079123427585245262,
        0.991455371120812639207),
      weights = c(
        0.0229353220105292249637,
        0.063092092629978553291,
        0.10479001032225018384,
        0.140653259715525918745,
        0.1690047266392679028266,
        0.1903505780647854099133,
        0.204432940075298892414,
        0.209482141084727828013,
        0.204432940075298892414,
        0.1903505780647854099133,
        0.169004726639267902827,
        0.140653259715525918745,
        0.1047900103222501838399,
        0.063092092629978553291,
        0.0229353220105292249637))
  } else if (nodes == 11) {
    list(
      points = c(
        -0.984085360094842464496,
        -0.906179845938663992798,
        -0.754166726570849220441,
        -0.5384693101056830910363,
        -0.2796304131617831934135,
        0,
        0.2796304131617831934135,
        0.5384693101056830910363,
        0.754166726570849220441,
        0.906179845938663992798,
        0.984085360094842464496),
      weights = c(
        0.042582036751081832865,
        0.1152333166224733940246,
        0.186800796556492657468,
        0.2410403392286475866999,
        0.272849801912558922341,
        0.2829874178574912132043,
        0.272849801912558922341,
        0.241040339228647586701,
        0.186800796556492657467,
        0.115233316622473394025,
        0.042582036751081832865))
  } else if (nodes == 7) {
    list(
      points = c(
        -0.9604912687080202834235,
        -0.7745966692414833770359,
        -0.4342437493468025580021,
        0,
        0.4342437493468025580021,
        0.7745966692414833770359,
        0.9604912687080202834235),
      weights = c(
        0.1046562260264672651938,
        0.268488089868333440729,
        0.401397414775962222905,
        0.450916538658474142345,
        0.401397414775962222905,
        0.268488089868333440729,
        0.104656226026467265194))
  } else stop("'qnodes' must be either 7, 11 or 15.")
}
