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
  eids <- ids  [event] # subject ids

  # quadrature points & weights, evaluated for each row of data
  qpts <- uapply(qp, unstandardise_qpts, t_beg, t_end)
  qwts <- uapply(qw, unstandardise_qwts, t_beg, t_end)
  qids <- rep(ids, qnodes)

  # quadrature points & weights, evaluated at upper limit of rows w/ interval censoring
  if (nicens) {
    ipts <- uapply(qp, unstandardise_qpts, t_beg[icens], t_upper)
    iwts <- uapply(qw, unstandardise_qwts, t_beg[icens], t_upper)
    iids <- rep(ids[icens], qnodes)
  } else {
    ipts <- rep(0,0)
    iwts <- rep(0,0)
    iids <- rep(0,0)
  }

  cpts <- c(epts, qpts, ipts)
  cids <- c(eids, qids, iids)

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
  t_state = data[ ,time_start]

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
        basis_epts,
        basis_qpts,
        basis_ipts,
        x,
        x_cpts,
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

  # get time start for initial, intermidiate and absorving states
  t_state = data[ ,time_start]

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
        t_state,
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


# Construct a list with information about the baseline hazard
#
# @param basehaz A string specifying the type of baseline hazard
# @param basehaz_ops A named list with elements df, knots
# @param ok_basehaz A list of admissible baseline hazards
# @param times A numeric vector with eventtimes for each individual
# @param status A numeric vector with event indicators for each individual
# @param upper_times A numeric vector (or NULL) with the upper limit for any
#   observations with interval censoring.
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
handle_basehaz2 <- function(basehaz,
                           basehaz_ops,
                           ok_basehaz     = c("weibull", "bs", "piecewise"),
                           ok_basehaz_ops = c("df", "knots"),
                           times,
                           status,
                           upper_times) {

  if (!basehaz %in% ok_basehaz)
    stop2("'basehaz' should be one of: ", comma(ok_basehaz))

  if (!all(names(basehaz_ops) %in% ok_basehaz_ops))
    stop2("'basehaz_ops' can only include: ", comma(ok_basehaz_ops))

  tt <- times[(status == 1)] # uncensored event times

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

    bknots <- get_bknots(c(times, upper_times))
    iknots <- get_iknots(tt, df = df, iknots = knots)
    basis  <- get_basis(tt, iknots = iknots, bknots = bknots, type = "bs")
    nvars  <- ncol(basis)  # number of aux parameters, basis terms

  } else if (basehaz == "ms") {

    df    <- basehaz_ops$df
    knots <- basehaz_ops$knots

    if (!is.null(df) && !is.null(knots)) {
      stop2("Cannot specify both 'df' and 'knots' for the baseline hazard.")
    }

    if (is.null(df)) {
      df <- 5L # default df for B-splines, assuming no intercept
      # NB this is ignored if the user specified knots
    }

    bknots <- get_bknots(c(times, upper_times))
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

    bknots <- get_bknots(c(times, upper_times))
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


#--------------- Functions related to association structure

# Return a named list summarising the association for one submodel
#
# @param user_x A character vector or NULL, being the user input to the
#   assoc argument (for one submodel) in the stan_jm call
# @param y_mod A list returned by a call to handle_glmod
# @param id_var The name of the ID variable
# @param metastuff A list with meta information about the joint model.
# @return A list with information about the desired association structure.
handle_assoc <- function(user_assoc, user_lag, y_mod, meta) {

  M              <- meta$M
  id_var         <- meta$id_var
  has_icens      <- meta$has_icens
  ok_assoc       <- meta$ok_assoc
  ok_assoc_data  <- meta$ok_assoc_data
  ok_assoc_int   <- meta$ok_assoc_int
  ok_assoc_icens <- meta$ok_assoc_icens

  ok_inputs <- c(ok_assoc, paste0(ok_assoc_data, "_data"),
                 uapply(ok_assoc_int, paste0, "_", ok_assoc_int))

  disallowed <- c("muslope", "shared_b", "shared_coef")

  trimmed_assoc <- trim_assoc(user_assoc, ok_assoc_data, ok_assoc_int)

  if (not.null(user_assoc) && !all(trimmed_assoc %in% ok_inputs))
    stop2("An unsupported association type has been specified. The ",
          "'assoc' argument can only include the following association ",
          "types: ", comma(ok_assoc), ", as well as possible interactions ",
          "either between association terms or with observed data.")

  if (any(trimmed_assoc %in% disallowed))
    stop2("The following association structures are temporarily disallowed: ",
          "and may be reinstated in a future release: ", comma(disallowed))

  if (any(ulist(has_icens)) && !all(trimmed_assoc %in% ok_assoc_icens))
    stop2("When interval censoring is present, only the following association ",
          "structures are allowed:", comma(ok_assoc_icens))

  if (not.null(user_assoc) && !is.character(user_assoc))
    stop2("The 'assoc' argument should be a character vector or, for a ",
          "multivariate joint model, possibly a list of character vectors.")

  assoc <- sapply(ok_inputs, `%in%`, trimmed_assoc, simplify = FALSE)
  if (is.null(user_assoc)) {
    assoc$null <- TRUE
  } else {
    if (assoc$null && (length(user_assoc) > 1L))
      stop2("In assoc, 'null' cannot be specified in conjuction ",
            "with another association type.")
    STOP_combination_not_allowed(assoc, "etavalue", "muvalue")
    STOP_combination_not_allowed(assoc, "etaslope", "muslope")
    STOP_combination_not_allowed(assoc, "etaauc",   "muauc")
  }

  # Parse suffix specifying indices for shared random effects
  cnms <- y_mod$z$group_cnms
  cnms_id <- cnms[[id_var]] # names of random effect terms
  assoc$which_b_zindex <- parse_assoc_sharedRE("shared_b",    user_assoc,
                                               max_index = length(cnms_id), cnms_id)
  assoc$which_coef_zindex <- parse_assoc_sharedRE("shared_coef", user_assoc,
                                                  max_index = length(cnms_id), cnms_id)

  if (length(intersect(assoc$which_b_zindex, assoc$which_coef_zindex)))
    stop("The same random effects indices should not be specified in both ",
         "'shared_b' and 'shared_coef'. Specifying indices in 'shared_coef' ",
         "will include both the fixed and random components.", call. = FALSE)

  if (length(assoc$which_coef_zindex)) {
    if (length(cnms) > 1L)
      stop("'shared_coef' association structure cannot be used when there is ",
           "clustering at levels other than the individual-level.", call. = FALSE)
    b_nms <- names(assoc$which_coef_zindex)
    assoc$which_coef_xindex <- sapply(b_nms, function(y, beta_nms) {
      beta_match <- grep(y, beta_nms, fixed = TRUE)
      if (!length(beta_match)) {
        stop("In association structure 'shared_coef', no matching fixed effect ",
             "component could be found for the following random effect: ", y,
             ". Perhaps consider using 'shared_b' association structure instead.")
      } else if (length(beta_match) > 1L) {
        stop("Bug found: In association structure 'shared_coef', multiple ",
             "fixed effect components have been found to match the following ",
             "random effect: ", y)
      }
      beta_match
    }, beta_nms = colnames(y_mod$X$X))
  } else assoc$which_coef_xindex <- numeric(0)

  if (!identical(length(assoc$which_coef_zindex), length(assoc$which_coef_xindex)))
    stop("Bug found: the lengths of the fixed and random components of the ",
         "'shared_coef' association structure are not the same.")

  # Parse suffix specifying formula for interactions with data
  ok_inputs_data <- paste0(ok_assoc_data, "_data")
  assoc$which_formulas <- sapply(ok_inputs_data, parse_assoc_data, user_assoc, simplify = FALSE)

  # Parse suffix specifying indices for interactions between association terms
  ok_inputs_interactions <- unlist(lapply(ok_assoc_int, paste0, "_", ok_assoc_int))
  assoc$which_interactions <- sapply(ok_inputs_interactions, parse_assoc_interactions,
                                     user_assoc, max_index = M, simplify = FALSE)

  # Lag for association structure
  assoc$which_lag <- user_lag

  assoc
}

# Remove suffixes from the user inputted assoc argument
#
# @param x A character vector, being the user input to the
#   assoc argument in the stan_jm call
# @param ok_assoc_data A character vector specifying which types
#   of association terms are allowed to be interacted with data
# @param ok_assoc_int A character vector specifying which types
#   of association terms are allowed to be interacted with other
#   association terms
trim_assoc <- function(x, ok_assoc_data, ok_assoc_int) {
  x <- gsub("^shared_b\\(.*",    "shared_b",    x)
  x <- gsub("^shared_coef\\(.*", "shared_coef", x)
  for (i in ok_assoc_data)
    x <- gsub(paste0("^", i, "_data\\(.*"),    paste0(i, "_data"), x)
  for (i in ok_assoc_int) for (j in ok_assoc_int)
    x <- gsub(paste0("^", i, "_", j, "\\(.*"), paste0(i, "_", j),  x)
  x
}

# Parse the indices specified for shared random effects
#
# @param x A character string corresponding to one of the allowed
#   association structures for shared random effects
# @param user_x A character vector, being the user input to the assoc
#   argument in the stan_jm call
# @param max_index An integer specifying the total number of random effects
#   in the longitudinal submodel, and therefore the maximum allowed index for
#   the shared random effects
# @param cnms The names of the random effects corresponding to the
#   individual-level (id_var) of clustering
# @return A numeric vector specifying indices for the shared random effects
parse_assoc_sharedRE <- function(x, user_x, max_index, cnms) {
  val <- grep(paste0("^", x, ".*"), user_x, value = TRUE)
  if (length(val)) {
    val2 <- unlist(strsplit(val, x))[-1]
    if (length(val2)) {
      index <- tryCatch(eval(parse(text = paste0("c", val2))), error = function(e)
        stop("Incorrect specification of the '", x, "' association structure. ",
             "See Examples in help file.", call. = FALSE))
      if (any(index > max_index))
        stop(paste0("The indices specified for the '", x, "' association structure are ",
                    "greater than the number of subject-specific random effects."), call. = FALSE)
    } else index <- seq_len(max_index)
    names(index) <- cnms[index]
    return(index)
  } else numeric(0)
}

# Parse the formula for specifying a data interaction with an association term
#
# @param x A character string corresponding to one of the allowed
#   association structures for interactions with data, for example,
#   "etavalue_data" or "etaslope_data"
# @param user_x A character vector, being the user input to the assoc
#   argument in the stan_jm call
# @return The parsed formula (which can be used for constructing a
#   design matrix for interacting data with association type x) or NULL
parse_assoc_data <- function(x, user_x) {
  val <- grep(paste0("^", x, ".*"), user_x, value = TRUE)
  if (length(val)) {
    val2 <- unlist(strsplit(val, x))[-1]
    fm <- tryCatch(eval(parse(text = val2)), error = function(e)
      stop(paste0("Incorrect specification of the formula in the '", x,
                  "' association structure. See Examples in the help file."), call. = FALSE))
    if (!is(fm, "formula"))
      stop(paste0("Suffix to '", x, "' association structure should include ",
                  "a formula within parentheses."), call. = FALSE)
    if (identical(length(fm), 3L))
      stop(paste0("Formula specified for '", x, "' association structure should not ",
                  "include a response."), call. = FALSE)
    if (length(lme4::findbars(fm)))
      stop(paste0("Formula specified for '", x, "' association structure should only ",
                  "include fixed effects."), call. = FALSE)
    if (fm[[2L]] == 1)
      stop(paste0("Formula specified for '", x, "' association structure cannot ",
                  "be an intercept only."), call. = FALSE)
    return(fm)
  } else numeric(0)
}


# Parse the indices specified for interactions between association terms
#
# @param x A character string corresponding to one of the allowed
#   association structures
# @param user_x A character vector, being the user input to the assoc
#   argument in the stan_jm call
# @param max_index An integer specifying the maximum allowed index
# @return A numeric vector specifying indices
parse_assoc_interactions <- function(x, user_x, max_index) {
  val <- grep(paste0("^", x, ".*"), user_x, value = TRUE)
  if (length(val)) {
    val2 <- unlist(strsplit(val, x))[-1]
    if (length(val2)) {
      index <- tryCatch(eval(parse(text = paste0("c", val2))), error = function(e)
        stop("Incorrect specification of the '", x, "' association structure. It should ",
             "include a suffix with parentheses specifying the indices of the association ",
             "terms you want to include in the interaction. See Examples in the help file.", call. = FALSE))
      if (any(index > max_index))
        stop("The indices specified for the '", x, "' association structure ",
             "cannot be greater than the number of longitudinal submodels.", call. = FALSE)
      return(index)
    } else
      stop("Incorrect specification of the '", x, "' association structure. It should ",
           "include a suffix with parentheses specifying the indices of the association ",
           "terms you want to include in the interaction. See Examples in the help file.", call. = FALSE)
  } else numeric(0)
}

# Make sure that interactions between association terms (for example
# etavalue_etaslope or mu_value_muvalue etc) are always ordered so that
# the first listed association term is for the submodel with the smallest
# index. For example, etavalue1_etavalue2 NOT etavalue2_etavalue1. This
# is to ensure there is no replication such as including both
# etavalue1_etavalue2 AND etavalue2_etavalue1 when passing to Stan.
#
# @param assoc A two dimensional array with information about desired association
#   structure for the joint model (returned by a call to validate_assoc).
# @param ok_assoc_int A character vector, specifying which association
#   structures are allowed to be used in interactions
check_order_of_assoc_interactions <- function(assoc, ok_assoc_int) {
  M <- ncol(assoc)
  for (i in ok_assoc_int) {
    for (j in ok_assoc_int) {
      header <- paste0(i, "_", j)
      header_reversed <- paste0(j, "_", i)
      for (m in 1:M) {
        if (assoc[header,][[m]]) {
          indices <- assoc["which_interactions",][[m]][[header]]
          sel <- which(indices < m)
          if (length(sel)) {
            # Remove indices for submodels before the current submodel m
            new_indices <- indices[-sel]
            assoc["which_interactions", ][[m]][[header]] <- new_indices
            assoc[header,][[m]] <- (length(new_indices) > 0L)
            # Replace those indices by reversing the order of association terms
            for (k in indices[sel]) {
              assoc["which_interactions",][[k]][[header_reversed]] <-
                unique(c(assoc["which_interactions",][[k]][[header_reversed]], m))
              assoc[header_reversed,][[k]] <-
                (length(assoc["which_interactions",][[k]][[header_reversed]]) > 0L)
            }
          }
        }
      }
    }
  }
  assoc
}


# Get the information need for combining the information in lower-level units
# clustered within an individual, when the patient-level is not the only
# clustering level in the longitudinal submodel
#
# @param cnms The component names for a single longitudinal submodel
# @param flist The flist for a single longitudinal submodel
# @param id_var The name of the ID variable
# @param qnodes Integer specifying the number of qnodes being used for
#   the GK quadrature in the stan_jm call
# @param grp_assoc Character string specifying the association structure used
#   for combining information in the lower level units clustered within an
#   individual
# @return A named list with the following elements:
#   has_grp: logical specifying whether the submodel has a grouping factor
#     that is clustered with patients.
#   grp_var: the name of any grouping factor that is clustered with patients.
#   grp_assoc: the user input to the grp_assoc argument in the stan_jm call.
#   grp_freq: a named vector with the number of lower level units clustered
#     within each individual.
#   grp_list: a named list containing the unique names for the lower level
#     units clustered within each individual.
get_basic_grp_info <- function(y_mod, id_var) {
  cnms  <- y_mod[["z"]][["group_cnms"]]
  flist <- y_mod[["z"]][["group_list"]]
  cnms_nms <- names(cnms)
  tally <- xapply(cnms_nms, FUN = function(x)
    # within each ID, count the number of levels for the grouping factor x
    tapply(flist[[x]], flist[[id_var]], FUN = n_distinct))
  sel <- which(sapply(tally, function(x) !all(x == 1L)) == TRUE)
  has_grp <- as.logical(length(sel))
  if (!has_grp) {
    return(nlist(has_grp))
  } else {
    if (length(sel) > 1L)
      stop("There can only be one grouping factor clustered within 'id_var'.")
    grp_var <- cnms_nms[sel]
    return(nlist(has_grp, grp_var))
  }
}

get_extra_grp_info <- function(basic_info, flist, id_var, grp_assoc,
                               ok_assoc_grp = c("sum", "mean", "min", "max")) {
  has_grp <- basic_info$has_grp
  grp_var <- basic_info$grp_var
  if (!has_grp) { # no grouping factor clustered within patients
    return(basic_info)
  } else { # submodel has a grouping factor clustered within patients
    if (is.null(grp_var))
      stop2("Bug found: could not find 'grp_var' in basic_info.")
    if (is.null(grp_assoc))
      stop2("'grp_assoc' cannot be NULL when there is a grouping factor ",
            "clustered within patients.")
    if (!grp_assoc %in% ok_assoc_grp)
      stop2("'grp_assoc' must be one of: ", comma(ok_assoc_grp))

    # cluster and patient ids for each row of the z matrix
    factor_grp <- factor(flist[[grp_var]])
    factor_ids <- factor(flist[[id_var]])

    # num clusters within each patient
    grp_freq <- tapply(factor_grp, factor_ids, FUN = n_distinct, simplify = FALSE)
    grp_freq <- unlist(grp_freq)

    # unique cluster ids for each patient id
    grp_list <- tapply(factor_grp, factor_ids, FUN = unique, simplify = FALSE)

    basic_info <- nlist(has_grp, grp_var)
    extra_info <- nlist(grp_assoc, grp_freq, grp_list)
    return(c(basic_info, extra_info))
  }
}



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
handle_assocmod2 <- function(data, assoc, y_mod, e_mod, grp_stuff, meta, j) {

  if (!requireNamespace("data.table"))
    stop2("the 'data.table' package must be installed to use this function.")

  id_var   <- meta$id_var
  time_var <- meta$time_var

  # before turning data into a data.table (for a rolling merge against
  # the quadrature points) we want to make sure that the data does not
  # include any NAs for the predictors or assoc formula variables
  tt <- attr(y_mod$terms, "term.labels")
  tt <- c(tt, uapply(assoc[["which_formulas"]], all.vars))
  fm <- reformulate(tt, response = NULL)
  df <- get_all_vars(fm, data)
  df <- df[complete.cases(df), , drop = FALSE]
  df <- df[e_mod$assoc_obs, ]

  # declare df as a data.table for merging with quadrature points
  dt <- prepare_data_table(df,
                           id_var   = id_var,
                           time_var = time_var,
                           grp_var  = grp_stuff$grp_var) # grp_var may be NULL
  ids_ <- seq_along( promote_to_factor(unique(df[ ,id_var])) )

  time_state <-  e_mod$cpts + suppressMessages(dplyr::left_join(
    data.frame(id = e_mod$cids),
    data.frame(id = ids_,
               t = e_mod$t_state)
                              ))$t

  # design matrices for calculating association structure based on
  # (possibly lagged) eta, slope, auc and any interactions with data
  args <- list(use_function = make_assoc_parts_for_stan,
               newdata      = dt,
               y_mod        = y_mod,
               grp_stuff    = grp_stuff,
               meta_stuff   = meta,
               assoc        = assoc,
               ids          = e_mod$cids,
               times        = time_state,
               j            = j)
  return(args)

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
                             ids, times, j, ...) {

  if (!requireNamespace("data.table"))
    stop("the 'data.table' package must be installed to use this function")

  id_var     <- meta_stuff$id_var
  time_var   <- meta_stuff$time_var
  epsilon    <- meta_stuff$epsilon
  auc_qnodes <- meta_stuff$auc_qnodes[j]

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

# Return design matrices for the longitudinal submodel. This is
# designed to generate the design matrices evaluated at the GK
# quadrature points, because it uses a 'terms' object to generate
# the model frame, and that terms object should have been generated
# from the longitudinal submodel's model frame when it was evaluated
# at the observation times; i.e. the predvars and X_bar would have
# come from the design matrices at the observation times, not the
# quadrature points.
#
# @param newdata A data frame; the data for the longitudinal submodel
#   at the event and quadrature points.
# @param y_mod The list returned by handle_y_mod, containing info about
#   the longitudinal submodel evaluated at the observation (not quadrature)
#   times, for example, the x_bar means used for centering, the predvars
#   attribute for the longitudinal submodel formula, and so on.
# @param include_Zt Whether to include the sparse Zt matrix in the
#   returned parts.
make_assoc_parts_for_stan <- function(newdata, y_mod, include_Zt = TRUE) {

  # construct model frame using predvars
  formula     <- use_predvars(y_mod, keep_response = FALSE)
  data        <- as.data.frame(newdata)
  model_frame <- stats::model.frame(lme4::subbars(formula), data)

  # fe design matrices
  x_form <- lme4::nobars(formula)
  x <- model.matrix(x_form, model_frame)
  xtemp <- drop_intercept(x)
  x_bar <- y_mod$x$x_bar
  xtemp <- sweep(xtemp, 2, x_bar, FUN = "-")

  # re design matrices
  bars <- lme4::findbars(formula)
  if (length(bars) > 2L)
    stop2("A maximum of 2 grouping factors are allowed.")
  z_parts <- lapply(bars, split_at_bars)
  z_forms <- fetch(z_parts, "re_form")
  z <- lapply(z_forms, model.matrix, model_frame)
  group_vars <- fetch(z_parts, "group_var")
  group_list <- lapply(group_vars, function(x) factor(model_frame[[x]]))
  names(z) <- names(group_list) <- group_vars

  ret <- nlist(x, xtemp, z, group_list, group_vars) # return list

  # optionally add the sparse Zt matrix
  if (include_Zt)
    ret$Zt <- lme4::mkReTrms(bars, model_frame)$Zt

  ret
}


# Function to substitute variables in the formula of a fitted model
# with the corresponding predvars based on the terms object for the model.
# (This is useful since lme4::glFormula doesn't allow a terms object to be
# passed as the first argument instead of a model formula).
#
# @param mod A (g)lmer model object from which to extract the formula and terms
# @return A reformulated model formula with variables replaced by predvars
use_predvars <- function(mod, keep_response = TRUE) {
  fm <- formula(mod)
  ff <- lapply(attr(terms(mod, fixed.only  = TRUE), "variables"), deparse, 500)[-1]
  fr <- lapply(attr(terms(mod, random.only = TRUE), "variables"), deparse, 500)[-1]
  pf <- lapply(attr(terms(mod, fixed.only  = TRUE), "predvars"),  deparse, 500)[-1]
  pr <- lapply(attr(terms(mod, random.only = TRUE), "predvars"),  deparse, 500)[-1]
  if (!identical(c(ff, fr), c(pf, pr))) {
    for (j in 1:length(ff))
      fm <- gsub(ff[[j]], pf[[j]], fm, fixed = TRUE)
    for (j in 1:length(fr))
      fm <- gsub(fr[[j]], pr[[j]], fm, fixed = TRUE)
  }
  rhs <- fm[[length(fm)]]
  if (is(rhs, "call"))
    rhs <- deparse(rhs, 500L)
  if (keep_response && length(fm) == 3L) {
    fm <- reformulate(rhs, response = formula(mod)[[2L]])
  } else if (keep_response && length(fm) == 2L) {
    warning2("No response variable found, reformulating RHS only.")
    fm <- reformulate(rhs, response = NULL)
  } else {
    fm <- reformulate(rhs, response = NULL)
  }
  fm
}

# Function to calculate the number of association parameters in the model
#
# @param assoc A list of length M with information about the association structure
#   type for each submodel, returned by an mapply call to validate_assoc
# @param a_mod_stuff A list of length M with the design matrices related to
#   the longitudinal submodels in the GK quadrature, returned by an mapply
#   call to handle_assocmod
# @return Integer indicating the number of association parameters in the model
get_num_assoc_pars <- function(assoc, a_mod_stuff) {
  sel1 <- c("etavalue", "etaslope", "etaauc",
            "muvalue", "muslope", "muauc")
  sel2 <- c("which_b_zindex", "which_coef_zindex")
  sel3 <- c("which_interactions")
  K1 <- sum(as.integer(assoc[sel1,]))
  K2 <- length(unlist(assoc[sel2,]))
  K3 <- length(unlist(assoc[sel3,]))
  K4 <- sum(fetch_(a_mod_stuff, "K_data"))
  K1 + K2 + K3 + K4
}


# Return the list of pars for Stan to monitor
#
# @param standata The list of data to pass to Stan
# @param is_jm A logical
# @return A character vector
pars_to_monitor <- function(standata, is_jm = FALSE) {
  c(if (standata$M > 0 && standata$intercept_type[1]) "yAlpha1",
    if (standata$M > 1 && standata$intercept_type[2]) "yAlpha2",
    if (standata$M > 2 && standata$intercept_type[3]) "yAlpha3",
    if (standata$M > 0 && standata$yK[1]) "yBeta1",
    if (standata$M > 1 && standata$yK[2]) "yBeta2",
    if (standata$M > 2 && standata$yK[3]) "yBeta3",
    if (is_jm) "e_alpha",
    if (is_jm && standata$e_K) "e_beta",
    if (is_jm && standata$a_K) "a_beta",
    if (standata$bK1 > 0) "b1",
    if (standata$bK2 > 0) "b2",
    if (standata$M > 0 && standata$has_aux[1]) "yAux1",
    if (standata$M > 1 && standata$has_aux[2]) "yAux2",
    if (standata$M > 2 && standata$has_aux[3]) "yAux3",
    if (is_jm && length(standata$basehaz_nvars)) "e_aux",
    if (standata$prior_dist_for_cov == 2 && standata$bK1 > 0) "bCov1",
    if (standata$prior_dist_for_cov == 2 && standata$bK2 > 0) "bCov2",
    if (standata$prior_dist_for_cov == 1 && standata$len_theta_L) "theta_L",
    "mean_PPD")
}


#--------------- Functions related to generating initial values

# Create a function that can be used to generate the model-based initial values for Stan
#
# @param e_mod_stuff A list object returned by a call to the handle_coxmod function
# @param standata The data list that will be passed to Stan
get_prefit_inits2 <- function(init_fit, e_mod, prior_stuff, prior_aux_stuff, standata) {

  init_means <- rstan::get_posterior_mean(init_fit)
  init_nms   <- rownames(init_means)
  inits <- generate_init_function2(e_mod, prior_stuff, prior_aux_stuff)()

  sel_b1 <- grep(paste0("^z_bMat1\\."), init_nms)
  if (length(sel_b1))
    inits[["z_bMat1"]] <- matrix(init_means[sel_b1,], nrow = standata$bK1)

  sel_b2 <- grep(paste0("^z_bMat2\\."), init_nms)
  if (length(sel_b2))
    inits[["z_bMat2"]] <- matrix(init_means[sel_b2,], nrow = standata$bK2)

  sel_bC1 <- grep(paste0("^bCholesky1\\."), init_nms)
  if (length(sel_bC1) > 1) {
    inits[["bCholesky1"]] <- matrix(init_means[sel_bC1,], nrow = standata$bK1)
  } else if (length(sel_bC1) == 1) {
    inits[["bCholesky1"]] <- aa(init_means[sel_bC1,])
  }

  sel_bC2 <- grep(paste0("^bCholesky2\\."), init_nms)
  if (length(sel_bC2) > 1) {
    inits[["bCholesky2"]] <- matrix(init_means[sel_bC2,], nrow = standata$bK2)
  } else if (length(sel_bC1) == 1) {
    inits[["bCholesky2"]] <- aa(init_means[sel_bC2,])
  }

  sel <- c("yGamma1", "yGamma2", "yGamma3",
           "z_yBeta1", "z_yBeta2", "z_yBeta3",
           "yAux1_unscaled", "yAux2_unscaled", "yAux3_unscaled",
           "bSd1", "bSd2", "z_b", "z_T", "rho", "zeta", "tau",
           "yGlobal1", "yGlobal2", "yGlobal3",
           "yLocal1", "yLocal2", "yLocal3",
           "yMix1", "yMix2", "yMix3",
           "yOol1", "yOol2", "yOol3")
  for (i in sel) {
    sel_i <- grep(paste0("^", i, "\\."), init_nms)
    if (length(sel_i))
      inits[[i]] <- aa(init_means[sel_i,])
  }
  return(function() inits)
}

generate_init_function2 <- function(e_mod_stuff, prior_stuff, prior_aux_stuff) {

  # Initial values for intercepts, coefficients and aux parameters
  if (e_mod_stuff$surv_type %in% c("right", "counting")) {
    e_beta <- e_mod_stuff$mod$coef
  } else if (e_mod_stuff$surv_type %in% c("interval", "interval2")) {
    e_beta <- -drop_intercept(e_mod_stuff$mod$coef) * e_mod_stuff$mod$scale
  } else {
    stop("Bug found: Invalid Surv type.")
  }
  e_aux <- if (e_mod_stuff$basehaz$type == 1L) runif(1, 0.5, 3) else rep(0, e_mod_stuff$basehaz$nvars)
  e_z_beta       <- standardise_coef(x        = e_beta,
                                     location = prior_stuff$prior_mean,
                                     scale    = prior_stuff$e_prior_scale)
  e_aux_unscaled <- standardise_coef(x        = e_aux,
                                     location = prior_aux_stuff$prior_mean_for_aux,
                                     scale    = prior_aux_stuff$prior_scale_for_aux)

  # Function to generate model based initial values
  model_based_inits <- rm_null(list(
    e_z_beta       = array_else_double(e_z_beta),
    e_aux_unscaled = array_else_double(e_aux_unscaled),
    e_gamma        = array_else_double(rep(0, e_mod_stuff$has_intercept))))

  return(function() model_based_inits)
}


# Function to construct a design matrix for the association structure in
# the event submodel, to be multiplied by a vector of association parameters
#
# @param assoc An array with information about the desired association
#   structure, returned by a call to validate_assoc.
# @param parts A list equal in length to the number of markers. Each element
#   parts[[m]] should contain a named list with components $mod_eta, $mod_eps,
#   $mod_auc, etc, which each contain either the linear predictor at quadtimes,
#   quadtimes + eps, and auc quadtimes, or the design matrices
#   used for constructing the linear predictor. Each element parts[[m]] should
#   also contain $X_data and $K_data.
# @param family A list of family objects, equal in length to the number of
#   longitudinal submodels.
# @param ... If parts does not contain the linear predictors, then this should
#   include elements beta and b, each being a length M list of parameters for the
#   longitudinal submodels.
# @return A design matrix containing the association terms to be multiplied by
#   the association paramters.
make_assoc_terms <- function(parts, assoc, family, ...) {
  M <- length(parts)
  a_X <- list()
  mark <- 1
  for (m in 1:M) {
    times   <- attr(parts[[m]], "times")
    epsilon <- attr(parts[[m]], "epsilon")
    qnodes  <- attr(parts[[m]], "auc_qnodes")
    qwts    <- attr(parts[[m]], "auc_qwts")

    eps_uses_derivative_of_x <-
      attr(parts[[m]], "eps_uses_derivative_of_x") # experimental

    has_assoc <- !assoc["null",][[m]]

    if (has_assoc) {
      assoc_m   <- assoc[,m]
      invlink_m <- family[[m]]$linkinv
      eta_m    <- get_element(parts, m = m, "eta", ...)
      eps_m    <- get_element(parts, m = m, "eps", ...)
      auc_m    <- get_element(parts, m = m, "auc", ...)
      X_data_m <- get_element(parts, m = m, "X_data", ...)
      K_data_m <- get_element(parts, m = m, "K_data", ...)
      grp_m    <- get_element(parts, m = m, "grp_stuff", ...)

      has_grp   <- grp_m$has_grp # TRUE/FALSE
      if (has_grp) {
        # method for collapsing information across clusters within patients
        grp_assoc <- grp_m$grp_assoc
        # indexing for collapsing across grps (based on the ids and times
        # used to generate the design matrices in make_assoc_parts)
        grp_idx <- attr(parts[[m]], "grp_idx")
      }

      #---  etavalue and any interactions  ---#

      # etavalue
      if (assoc_m[["etavalue"]]) {
        if (has_grp) {
          a_X[[mark]] <- collapse_within_groups(eta_m, grp_idx, grp_assoc)
        } else {
          a_X[[mark]] <- eta_m
        }
        mark <- mark + 1
      }

      # etavalue * data interactions
      if (assoc_m[["etavalue_data"]]) {
        X_temp <- X_data_m[["etavalue_data"]]
        K_temp <- K_data_m[["etavalue_data"]]
        for (i in 1:K_temp) {
          if (is.matrix(eta_m)) {
            val <- sweep(eta_m, 2L, X_temp[, i], `*`)
          } else {
            val <- as.vector(eta_m) * X_temp[, i]
          }
          if (has_grp) {
            a_X[[mark]] <- collapse_within_groups(val, grp_idx, grp_assoc)
          } else {
            a_X[[mark]] <- val
          }
          mark <- mark + 1
        }
      }

      # etavalue * etavalue interactions
      if (assoc_m[["etavalue_etavalue"]]) {
        sel <- assoc_m[["which_interactions"]][["etavalue_etavalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          val <- eta_m * eta_j
          a_X[[mark]] <- val
          mark <- mark + 1
        }
      }

      # etavalue * muvalue interactions
      if (assoc_m[["etavalue_muvalue"]]) {
        sel <- assoc_m[["which_interactions"]][["etavalue_muvalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          invlink_j <- family[[j]]$linkinv
          val <- eta_m * invlink_j(eta_j)
          a_X[[mark]] <- val
          mark <- mark + 1
        }
      }

      #---  etaslope and any interactions  ---#

      if (assoc_m[["etaslope"]] || assoc_m[["etaslope_data"]]) {
        if (eps_uses_derivative_of_x) {
          deta_m <- eps_m
        } else {
          deta_m <- (eps_m - eta_m) / epsilon
        }
      }

      # etaslope
      if (assoc_m[["etaslope"]]) {
        if (has_grp) {
          a_X[[mark]] <- collapse_within_groups(deta_m, grp_idx, grp_assoc)
        } else {
          a_X[[mark]] <- deta_m
        }
        mark <- mark + 1
      }

      # etaslope * data interactions
      if (assoc_m[["etaslope_data"]]) {
        X_temp <- X_data_m[["etaslope_data"]]
        K_temp <- K_data_m[["etaslope_data"]]
        for (i in 1:K_temp) {
          if (is.matrix(deta_m)) {
            val <- sweep(deta_m, 2L, X_temp[, i], `*`)
          } else {
            val <- as.vector(deta_m) * X_temp[, i]
          }
          if (has_grp) {
            a_X[[mark]] <- collapse_within_groups(val, grp_idx, grp_assoc)
          } else {
            a_X[[mark]] <- val
          }
          mark <- mark + 1
        }
      }

      #---  etaauc  ---#

      if (assoc_m[["etaauc"]]) {
        if (is.matrix(eta_m)) {
          nr <- nrow(eta_m)
          nc <- ncol(eta_m)
          val   <- matrix(NA, nrow = nr, ncol = nc)
          for (j in 1:nc) {
            wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
            auc_j <- auc_m[, ((j-1) * qnodes + 1):(j * qnodes), drop = FALSE]
            tmp_j <- sweep(auc_j, 2L, wgt_j, `*`)
            val[,j] <- rowSums(tmp_j)
          }
        } else {
          val <- c()
          for (j in 1:length(eta_m)) {
            wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
            auc_j <- auc_m[((j-1) * qnodes + 1):(j * qnodes)]
            val[j] <- sum(wgt_j * auc_j)
          }
        }
        a_X[[mark]] <- val
        mark <- mark + 1
      }

      #---  muvalue and any interactions  ---#

      # muvalue
      if (assoc_m[["muvalue"]]) {
        mu_m <- invlink_m(eta_m)
        a_X[[mark]] <- mu_m
        mark <- mark + 1
      }

      # muvalue * data interactions
      if (assoc_m[["muvalue_data"]]) {
        mu_m <- invlink_m(eta_m)
        X_temp <- X_data_m[["muvalue_data"]]
        K_temp <- K_data_m[["muvalue_data"]]
        for (i in 1:K_temp) {
          if (is.matrix(mu_m)) {
            val <- sweep(mu_m, 2L, X_temp[, i], `*`)
          } else {
            val <- as.vector(mu_m) * X_temp[, i]
          }
          if (has_grp) {
            a_X[[mark]] <- collapse_within_groups(val, grp_idx, grp_assoc)
          } else {
            a_X[[mark]] <- val
          }
          mark <- mark + 1
        }
      }

      # muvalue * etavalue interactions
      if (assoc_m[["muvalue_etavalue"]]) {
        sel <- assoc_m[["which_interactions"]][["muvalue_etavalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          val   <- invlink_m(eta_m) * eta_j
          a_X[[mark]] <- val
          mark <- mark + 1
        }
      }

      # muvalue * muvalue interactions
      if (assoc_m[["muvalue_muvalue"]]) {
        sel <- assoc_m[["which_interactions"]][["muvalue_muvalue"]]
        for (j in sel) {
          eta_j <- get_element(parts, m = j, "eta", ...)
          invlink_j <- family[[j]]$linkinv
          val <- invlink_m(eta_m) * invlink_j(eta_j)
          a_X[[mark]] <- val
          mark <- mark + 1
        }
      }

      #---  muslope and any interactions  ---#

      if (assoc_m[["muslope"]] || assoc_m[["muslope_data"]]) {
        if (eps_uses_derivative_of_x) {
          stop2("Cannot currently use muslope interaction structure.")
        } else {
          dmu_m <- (invlink_m(eps_m) - invlink_m(eta_m)) / epsilon
        }
      }

      # muslope
      if (assoc_m[["muslope"]]) {
        a_X[[mark]] <- dmu_m
        mark <- mark + 1
      }

      # muslope * data interactions
      if (assoc_m[["muslope_data"]]) {
        X_temp <- X_data_m[["muslope_data"]]
        K_temp <- K_data_m[["muslope_data"]]
        for (i in 1:K_temp) {
          if (is.matrix(dmu_m)) {
            val <- sweep(dmu_m, 2L, X_temp[, i], `*`)
          } else {
            val <- as.vector(dmu_m) * X_temp[, i]
          }
          if (has_grp) {
            a_X[[mark]] <- collapse_within_groups(val, grp_idx, grp_assoc)
          } else {
            a_X[[mark]] <- val
          }
          mark <- mark + 1
        }
      }

      #---  muauc  ---#

      if (assoc_m[["muauc"]]) {
        if (is.matrix(eta_m)) {
          nr <- nrow(eta_m)
          nc <- ncol(eta_m)
          val   <- matrix(NA, nrow = nr, ncol = nc)
          for (j in 1:nc) {
            wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
            auc_j <- invlink_m(auc_m[, ((j-1) * qnodes + 1):(j * qnodes), drop = FALSE])
            tmp_j <- sweep(auc_j, 2L, wgt_j, `*`)
            val[,j] <- rowSums(tmp_j)
          }
        } else {
          val <- c()
          for (j in 1:length(eta_m)) {
            wgt_j <- qwts[((j-1) * qnodes + 1):(j * qnodes)]
            auc_j <- invlink_m(auc_m[((j-1) * qnodes + 1):(j * qnodes)])
            val[j] <- sum(wgt_j * auc_j)
          }
        }
        a_X[[mark]] <- val
        mark <- mark + 1
      }

    }
  }
  for (m in 1:M) {
    # shared_b
    if (assoc["shared_b",][[m]]) {
      sel <- assoc["which_b_zindex",][[m]]
      val <- get_element(parts, m = m, "b_mat", ...)[,sel]
      a_X[[mark]] <- val
      mark <- mark + 1
    }
  }
  for (m in 1:M) {
    # shared_coef
    if (assoc["shared_coef",][[m]]) {
      sel <- assoc["which_coef_zindex",][[m]]
      val <- get_element(parts, m = m, "b_mat", ...)[,sel]
      a_X[[mark]] <- val
      mark <- mark + 1
    }
  }

  if (is.matrix(a_X[[1L]])) a_X else do.call("cbind", a_X)
}

# Function to get an "element" (e.g. a linear predictor, a linear predictor
# evaluated at epsilon shift, linear predictor evaluated at auc quadpoints,
# etc) constructed from the "parts" (e.g. mod_eta, mod_eps, mod_auc, etc)
# returned by a call to the function 'make_assoc_parts'.
#
# @param parts A named list containing the parts for constructing the association
#   structure. It may contain elements $mod_eta, $mod_eps, $mod_auc, etc. as
#   well as $X_data, $K_data, $grp_stuff. It is returned by a call to the
#   function 'make_assoc_parts'.
# @param m An integer specifying which submodel to get the element for.
# @param which A character string specifying which element to get.
get_element <- function(parts, m = 1, which = "eta", ...) {

  ok_which_args <- c("eta", "eps", "auc", "X_data", "K_data",
                     "b_mat", "grp_stuff")
  if (!which %in% ok_which_args)
    stop("'which' must be one of: ", paste(ok_which_args, collapse = ", "))

  if (which %in% c("eta", "eps", "auc")) {
    part <- parts[[m]][[paste0("mod_", which)]]
    if (is.null(part)) {
      # model doesn't include an assoc related to 'which'
      return(NULL)
    } else {
      # construct linear predictor for the 'which' part
      x <- part$x
      Zt <- part$Zt
      Znames  <- part$Z_names
      if (is.null(x) || is.null(Zt))
        stop2("Bug found: cannot find x and Zt in 'parts'. They are ",
              "required to build the linear predictor for '", which, "'.")

      dots <- list(...)
      beta <- dots$beta[[m]]
      b    <- dots$b[[m]]
      if (is.null(beta) || is.null(b))
        stop2("Bug found: beta and b must be provided to build the ",
              "linear predictor for '", which, "'.")

      eta <- linear_predictor(beta, x)
      if (NCOL(b) == 1) {
        eta <- eta + as.vector(b %*% Zt)
      } else {
        eta <- eta + as.matrix(b %*% Zt)
      }
      return(eta)
    }
  } else if (which %in% c("X_data", "K_data", "b_mat", "grp_stuff")) {
    return(parts[[m]][[which]])
  } else {
    stop("'which' argument doesn't include a valid entry.")
  }
}

# Collapse the linear predictor across the lower level units
# clustered an individual, using the function specified in the
# 'grp_assoc' argument
#
# @param eta The linear predictor evaluated for all lower level groups
#   at the quadrature points.
# @param grp_idx An N*2 array providing the indices of the first (col 1)
#   and last (col 2) observations in eta that correspond to individuals
#   i = 1,...,N.
# @param grp_assoc Character string, the function to use to collapse
#   across the lower level units clustered within individuals.
# @return A vector or matrix, depending on the method called.
collapse_within_groups <- function(eta, grp_idx, grp_assoc = "sum") {
  UseMethod("collapse_within_groups")
}
collapse_within_groups.default <- function(eta, grp_idx, grp_assoc) {
  N <- nrow(grp_idx)
  val <- rep(NA, N)
  for (n in 1:N) {
    tmp <- eta[grp_idx[n,1]:grp_idx[n,2]]
    val[n] <- do.call(grp_assoc, list(tmp))
  }
  val
}
collapse_within_groups.matrix <- function(eta, grp_idx, grp_assoc) {
  N <- nrow(grp_idx)
  val <- matrix(NA, nrow = nrow(eta), ncol = N)
  for (n in 1:N) {
    tmp <- eta[, grp_idx[n,1]:grp_idx[n,2], drop = FALSE]
    val[,n] = apply(tmp, 1L, grp_assoc)
  }
  val
}
