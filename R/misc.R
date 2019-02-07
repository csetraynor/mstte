#------------ General
# Check whether a vector/matrix/array contains an "(Intercept)"
check_for_intercept <- function(x, logical = FALSE) {
  nms <- if (is.matrix(x)) colnames(x) else names(x)
  sel <- which("(Intercept)" %in% nms)
  if (logical) as.logical(length(sel)) else sel
}

# Drop intercept from a vector/matrix/array of named coefficients
drop_intercept <- function(x) {
  sel <- check_for_intercept(x)
  if (length(sel) && is.matrix(x)) {
    x[, -sel, drop = FALSE]
  } else if (length(sel)) {
    x[-sel]
  } else {
    x
  }
}

# shorthand for as array and unlist
aau <- function(x){
  aa(ulist(x))
}

# shorthand for do.call rbind
rb <- function(x){
  do.call(rbind, x)
}

# Select rows of a matrix
#
# @param x A matrix.
# @param rows Logical or numeric vector stating which rows of 'x' to retain.
keep_rows <- function(x, rows = 1:nrow(x)) {
  x[rows, , drop = FALSE]
}

# Concatenate (i.e. 'c(...)') but don't demote factors to integers
ulist <- function(...) { unlist(list(...)) }

# Create a named list using specified names or, if names are omitted, using the
# names of the objects in the list
#
# @param ... Objects to include in the list.
# @return A named list.
nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name))
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  }

  return(out)
}

# Unlist the result from an lapply call
#
# @param X,FUN,... Same as lapply
uapply <- function(X, FUN, ...) {
  unlist(lapply(X, FUN, ...))
}

# Unlist the result from an lapply call not recursive
#
# @param X,FUN,... Same as lapply
nonrecuapply <- function(X, FUN, ...) {
  unlist(lapply(X, FUN, ...), recursive = FALSE)
}

# Stop without printing call
stop2    <- function(...) stop(..., call. = FALSE)

# Immediate warning without printing call
warning2 <- function(...) warning(..., immediate. = TRUE, call. = FALSE)

# Shorthand for suppress warnings
SW <- function(expr) base::suppressWarnings(expr)

# Check if an object is NULL
is_null <- function(x) {
  is.null(x) || ifelse(is.vector(x), all(sapply(x, is.null)), FALSE)
}

# Recursively removes NULL entries from an object
rm_null <- function(x, recursive = TRUE) {
  x <- Filter(Negate(is_null), x)
  if (recursive) {
    x <- lapply(x, function(x) if (is.list(x)) rm_null(x) else x)
  }
  x
}

# Shorthand for as.integer, as.double, as.matrix, as.array
ai <- function(x, ...) as.integer(x, ...)
ad <- function(x, ...) as.double (x, ...)
am <- function(x, ...) as.matrix (x, ...)
aa <- function(x, ...) as.array  (x, ...)
# Safe deparse
safe_deparse <- function(expr) deparse(expr, 500L)

# Maybe broadcast
#
# @param x A vector or scalar.
# @param n Number of replications to possibly make.
# @return If \code{x} has no length the \code{0} replicated \code{n} times is
#   returned. If \code{x} has length 1, the \code{x} replicated \code{n} times
#   is returned. Otherwise \code{x} itself is returned.
maybe_broadcast <- function(x, n) {
  if (!length(x)) {
    rep(0, times = n)
  } else if (length(x) == 1L) {
    rep(x, times = n)
  } else {
    x
  }
}

# Return the cutpoints for a specified number of quantiles of 'x'
#
# @param x A numeric vector.
# @param nq Integer specifying the number of quantiles.
# @return A vector of percentiles corresponding to percentages 100*k/m for
#   k=1,2,...,nq-1.
qtile <- function(x, nq = 2) {
  if (nq > 1) {
    probs <- seq(1, nq - 1) / nq
    return(quantile(x, probs = probs))
  } else if (nq == 1) {
    return(NULL)
  } else {
    stop("'nq' must be >= 1.")
  }
}

# -------------- Sanity checks ---------#

# Issue warning if high rhat values
#
# @param rhats Vector of rhat values.
# @param threshold Threshold value. If any rhat values are above threshold a
#   warning is issued.
check_rhats <- function(rhats, threshold = 1.1, check_lp = FALSE) {
  if (!check_lp)
    rhats <- rhats[!names(rhats) %in% c("lp__", "log-posterior")]

  if (any(rhats > threshold, na.rm = TRUE))
    warning("Markov chains did not converge! Do not analyze results!",
            call. = FALSE, noBreaks. = TRUE)
}

# Check that a tmat is concordant with
# a formula list
check_trans <- function(f, t){
  if(length(f) != sum(!is.na(t)) ){
    stop2("Formula list has to include a formula for each transition as expressed in the transition matrix")
  }
}

# Throw error if parameter isn't a positive scalar
#
# @param x The object to test.
validate_positive_scalar <- function(x, not_greater_than = NULL) {
  nm <- deparse(substitute(x))
  if (is.null(x))
    stop(nm, " cannot be NULL", call. = FALSE)
  if (!is.numeric(x))
    stop(nm, " should be numeric", call. = FALSE)
  if (any(x <= 0))
    stop(nm, " should be postive", call. = FALSE)
  if (!is.null(not_greater_than)) {
    if (!is.numeric(not_greater_than) || (not_greater_than <= 0))
      stop("'not_greater_than' should be numeric and postive")
    if (!all(x <= not_greater_than))
      stop(nm, " should less than or equal to ", not_greater_than, call. = FALSE)
  }
}

# Check if priors were autoscaled
#
# @param prior_stuff A list with prior info returned by handle_glm_prior
# @param has A logical checking, for example, whether the model has_predictors,
#   has_intercept, has_assoc, etc
# @param adjusted_prior_scale The prior scale after any autoscaling
check_if_rescaled <- function(prior_stuff, has, adjusted_prior_scale) {
  prior_stuff$prior_autoscale && has &&
    !is.na(prior_stuff$prior_dist_name) &&
    !all(prior_stuff$prior_scale == adjusted_prior_scale)
}

# Check that a stanfit object (or list returned by rstan::optimizing) is valid
#
check_stanfit <- function(x) {
  if (is.list(x)) {
    if (!all(c("par", "value") %in% names(x)))
      stop("Invalid object produced please report bug")
  }
  else {
    stopifnot(is(x, "stanfit"))
    if (x@mode != 0)
      stop("Invalid stanfit object produced please report bug")
  }
  return(TRUE)
}

# Test if an object inherits a specific stanreg class
#
# @param x The object to test.
is.stanreg   <- function(x) inherits(x, "stanreg")
is.stansurv  <- function(x) inherits(x, "stansurv")
is.stanmvreg <- function(x) inherits(x, "stanmvreg")
is.stanjm    <- function(x) inherits(x, "stanjm")
is.stanmstte <- function(x) inherits(x, "stanmstte")

# ------------- Helpers ---------------#
#Compute point estimates and standard errors from pointwise vectors
#
# @param x A matrix.
# @return An ncol(x) by 2 matrix with columns 'Estimate' and 'SE'
#   and rownames equal to colnames(x).
#
table_of_estimates <- function(x) {
  out <- cbind(
    Estimate = colSums2(x),
    SE = sqrt(nrow(x) * colVars(x))
  )
  rownames(out) <- colnames(x)
  return(out)
}

append_trans <- function(x, i, labs){
  if(is.null(labs)){
    paste0(x," trans(",i,")" )
  } else {
    paste0(x," trans(",labs,")" )
  }
}

append_title <- function(object){
  if(is.null(object$transition_labels)){
    paste("\nTransition", seq_len(object$length_of_transition),"\n")
  } else {
    paste("\nTransition", object$transition_labels,"\n")
  }
}

get_transition_name <- function(obj, i){

  append_trans(NULL, i, obj$transition_labels[i])
}

# Return the name for the intercept parameter
get_int_name_basehaz <- function(x, is_jm = FALSE, ...) {
  if (is_jm || has_intercept(x)) "(Intercept)" else NULL
}

# Return the names for the auxiliary parameters
get_aux_name_basehaz <- function(x, ...) {
  switch(get_basehaz_name(x),
         exp       = NULL,
         weibull   = "weibull-shape",
         gompertz  = "gompertz-scale",
         ms        = paste0("m-splines-coef", seq(x$nvars)),
         bs        = paste0("b-splines-coef", seq(x$nvars)),
         piecewise = paste0("piecewise-coef", seq(x$nvars)),
         NA)
}

#' Extract X, Y or Z from a stanreg object
#'
#' @keywords internal
#' @export
#' @param ... Other arguments passed to methods. For a \code{stanmvreg} object
#'   this can be an integer \code{m} specifying the submodel.
#' @return For \code{get_x} and \code{get_z}, a matrix. For \code{get_y}, either
#'   a vector or a matrix, depending on how the response variable was specified.
get_y <- function(object, ...) UseMethod("get_y")
#' @rdname get_y
#' @export
get_x <- function(object, ...) UseMethod("get_x")
#' @export
get_y.default <- function(object, ...) {
  object[["y"]] %ORifNULL% model.response(model.frame(object))
}
#' @export
get_x.default <- function(object, ...) {
  object[["x"]] %ORifNULL% model.matrix(object)
}
#' @export
get_y.stanmstte <- function(object, ind, ...) {
  object[["y"]][[ind]] %ORifNULL% model.response(model.frame(object))
}
#' @export
get_x.default <- function(object, ind, ...) {
  object[["x"]][[ind]] %ORifNULL% model.matrix(object)
}

# Wrapper for rstan::summary
# @param stanfit A stanfit object created using rstan::sampling or rstan::vb
# @return A matrix of summary stats
make_stan_summary <- function(stanfit) {
  levs <- c(0.5, 0.8, 0.95)
  qq <- (1 - levs) / 2
  probs <- sort(c(0.5, qq, 1 - qq))
  rstan::summary(stanfit, probs = probs, digits = 10)$summary
}

# Get the correct column name to use for selecting the median
#
# @param algorithm String naming the estimation algorithm (probably
#   \code{fit$algorithm}).
# @return Either \code{"50%"} or \code{"Median"} depending on \code{algorithm}.
select_median <- function(algorithm) {
  switch(algorithm,
         sampling = "50%",
         meanfield = "50%",
         fullrank = "50%",
         optimizing = "Median",
         stop("Bug found (incorrect algorithm name passed to select_median)",
              call. = FALSE))
}
# Replace the parameter names slot of an object of class 'stanfit'.
#
# @param stanfit An object of class 'stanfit'.
# @param new_nms A character vector of new parameter names.
# @return A 'stanfit' object.
replace_stanfit_nms <- function(stanfit, new_nms) {
  stanfit@sim$fnames_oi <- new_nms
  stanfit
}

# Check formula object
#
# @param formula The user input to the formula argument.
# @param needs_response A logical; if TRUE then formula must contain a LHS.
validate_formula <- function(formula, needs_response = TRUE) {

  if (!inherits(formula, "formula")) {
    stop2("'formula' must be a formula.")
  }

  if (needs_response) {
    len <- length(formula)
    if (len < 3) {
      stop2("'formula' must contain a response.")
    }
  }
  as.formula(formula)
}

# Check object is a Surv object with a valid type
#
# @param x A Surv object; the LHS of a formula evaluated in a data frame environment.
# @param ok_types A character vector giving the allowed types of Surv object.
validate_surv <- function(x, ok_types = c("right", "counting",
                                          "interval", "interval2")) {
  if (!inherits(x, "Surv"))
    stop2("LHS of 'formula' must be a 'Surv' object.")
  if (!attr(x, "type") %in% ok_types)
    stop2("Surv object type must be one of: ", comma(ok_types))
  x
}

# If a is NULL (and Inf, respectively) return b, otherwise just return a
# @param a,b Objects
`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}
`%ORifINF%` <- function(a, b) {
  if (a == Inf) b else a
}


# Return the list with summary information about the baseline hazard
#
# @return A named list.
get_basehaz <- function(x, ind) {
  if (is.stansurv(x))
    return(x$basehaz)
  if (is.stanjm(x))
    return(x$survmod$basehaz)
  if (is.stanmstte(x))
    return(x$basehaz[[ind]])
  stop("Bug found: could not find basehaz.")
}



# ------- Helpers data ------------------ #
# Parse the model formula
#
# @param formula The user input to the formula argument.
# @param data The user input to the data argument (i.e. a data frame).
parse_formula <- function(formula, data) {

  formula <- validate_formula(formula, needs_response = TRUE)

  lhs      <- lhs(formula) # full LHS of formula
  lhs_form <- reformulate_lhs(lhs)

  rhs        <- rhs(formula)         # RHS as expression
  rhs_form   <- reformulate_rhs(rhs) # RHS as formula
  rhs_terms  <- terms(rhs_form, specials = "tde")
  rhs_vars   <- rownames(attr(rhs_terms, "factors"))

  allvars <- all.vars(formula)
  allvars_form <- reformulate(allvars)

  surv <- eval(lhs, envir = data) # Surv object
  surv <- validate_surv(surv)
  type <- attr(surv, "type")

  if (type == "right") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[3L]])
    min_t    <- 0
    max_t    <- max(surv[, "time"])
  } else if (type == "counting") {
    tvar_beg <- as.character(lhs[[2L]])
    tvar_end <- as.character(lhs[[3L]])
    dvar     <- as.character(lhs[[4L]])
    min_t    <- min(surv[, "start"])
    max_t    <- max(surv[, "stop"])
  } else if (type == "interval") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[4L]])
    min_t    <- 0
    max_t    <-  max(surv[, c("time1", "time2")])
  } else if (type == "interval2") {
    tvar_beg <- NULL
    tvar_end <- as.character(lhs[[2L]])
    dvar     <- as.character(lhs[[3L]])
    min_t    <- 0
    max_t    <- max(surv[, c("time1", "time2")])
  }

  sel <- attr(rhs_terms, "specials")$tde

  if (!is.null(sel)) { # model has tde

    # replace 'tde(x, ...)' in formula with 'x'
    tde_oldvars <- rhs_vars
    tde_newvars <- sapply(tde_oldvars, function(oldvar) {
      if (oldvar %in% rhs_vars[sel]) {
        tde <- function(newvar, ...) { # define tde function locally
          safe_deparse(substitute(newvar))
        }
        eval(parse(text = oldvar))
      } else oldvar
    }, USE.NAMES = FALSE)
    term_labels <- attr(rhs_terms, "term.labels")
    for (i in sel) {
      sel_terms <- which(attr(rhs_terms, "factors")[i, ] > 0)
      for (j in sel_terms) {
        term_labels[j] <- gsub(tde_oldvars[i],
                               tde_newvars[i],
                               term_labels[j],
                               fixed = TRUE)
      }
    }
    tf_form <- reformulate(term_labels, response = lhs)

    # extract 'tde(x, ...)' from formula and construct 'bs(times, ...)'
    tde_terms <- lapply(rhs_vars[sel], function(x) {
      tde <- function(vn, ...) { # define tde function locally
        dots <- list(...)
        ok_args <- c("df")
        if (!isTRUE(all(names(dots) %in% ok_args)))
          stop2("Invalid argument to 'tde' function. ",
                "Valid arguments are: ", comma(ok_args))
        df <- if (is.null(dots$df)) 3 else dots$df
        degree <- 3
        if (df == 3) {
          dots[["knots"]] <- numeric(0)
        } else {
          dx <- (max_t - min_t) / (df - degree + 1)
          dots[["knots"]] <- seq(min_t + dx, max_t - dx, dx)
        }
        dots[["Boundary.knots"]] <- c(min_t, max_t)
        sub("^list\\(", "bs\\(times__, ", safe_deparse(dots))
      }
      tde_calls <- eval(parse(text = x))
      sel_terms <- which(attr(rhs_terms, "factors")[x, ] > 0)
      new_calls <- sapply(seq_along(sel_terms), function(j) {
        paste0(term_labels[sel_terms[j]], ":", tde_calls)
      })
      nlist(tde_calls, new_calls)
    })
    td_basis <- fetch(tde_terms, "tde_calls")
    new_calls <- fetch_(tde_terms, "new_calls")
    td_form <- reformulate(new_calls, response = NULL, intercept = FALSE)

  } else { # model doesn't have tde
    tf_form  <- formula
    td_form  <- NULL
    td_basis <- NULL
  }

  nlist(formula,
        lhs,
        rhs,
        lhs_form,
        rhs_form,
        tf_form,
        td_form,
        td_basis,
        fe_form = rhs_form, # no re terms accommodated yet
        re_form = NULL,     # no re terms accommodated yet
        allvars,
        allvars_form,
        tvar_beg,
        tvar_end,
        dvar,
        surv_type = attr(surv, "type"))
}


# Reformulate as RHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
reformulate_rhs <- function(x) {
  #x <- deparse(x, 500L)
  x <- formula(substitute(~ RHS, list(RHS = x)))
  x
}

# Extract LHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
lhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[2L]]
  } else {
    out <- NULL
  }
  out
}

# Return the model frame
#
# @param formula The parsed model formula.
# @param data The model data frame.
make_model_frame <- function(formula, data, check_constant = TRUE) {

  # construct terms object from formula
  Terms <- terms(formula)

  # construct model frame
  mf <- model.frame(Terms, data)

  # check no constant vars
  if (check_constant)
    mf <- check_constant_vars(mf)

  # check for terms
  mt <- attr(mf, "terms")
  if (is.empty.model(mt))
    stop2("No intercept or predictors specified.")

  nlist(mf, mt)
}


# Return the response vector (time) for estimation
#
# @param model_frame The model frame.
# @param type The type of time variable to return:
#   "beg": the entry time for the row in the survival data,
#   "end": the exit  time for the row in the survival data,
#   "gap": the difference between entry and exit times,
#   "upp": if the row involved interval censoring, then the exit time
#          would have been the lower limit of the interval, and "upp"
#          is the upper limit of the interval.
# @return A numeric vector.
make_t <- function(model_frame, type = c("beg", "end", "gap", "upp")) {

  type <- match.arg(type)
  resp <- if (survival::is.Surv(model_frame))
    model_frame else model.response(model_frame)
  surv <- attr(resp, "type")
  err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")

  t_beg <- switch(surv,
                  "right"     = rep(0, nrow(model_frame)),
                  "interval"  = rep(0, nrow(model_frame)),
                  "interval2" = rep(0, nrow(model_frame)),
                  "counting"  = as.vector(resp[, "start"]),
                  stop(err))

  t_end <- switch(surv,
                  "right"     = as.vector(resp[, "time"]),
                  "interval"  = as.vector(resp[, "time1"]),
                  "interval2" = as.vector(resp[, "time1"]),
                  "counting"  = as.vector(resp[, "stop"]),
                  stop(err))

  t_upp <- switch(surv,
                  "right"     = rep(NaN, nrow(model_frame)),
                  "counting"  = rep(NaN, nrow(model_frame)),
                  "interval"  = as.vector(resp[, "time2"]),
                  "interval2" = as.vector(resp[, "time2"]),
                  stop(err))

  switch(type,
         "beg" = t_beg,
         "end" = t_end,
         "gap" = t_end - t_beg,
         "upp" = t_upp,
         stop("Bug found: cannot handle specified 'type'."))
}


# Return the response vector (status indicator)
#
# @param model_frame The model frame.
# @return A numeric vector.
make_d <- function(model_frame) {

  resp <- if (survival::is.Surv(model_frame))
    model_frame else model.response(model_frame)
  surv <- attr(resp, "type")
  err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")

  switch(surv,
         "right"     = as.vector(resp[, "status"]),
         "interval"  = as.vector(resp[, "status"]),
         "interval2" = as.vector(resp[, "status"]),
         "counting"  = as.vector(resp[, "status"]),
         stop(err))
}

# Check if any variables in a model frame are constants
# (the exception is that a constant variable of all 1's is allowed)
#
# @param mf A model frame or model matrix
# @return If no constant variables are found mf is returned, otherwise an error
#   is thrown.
check_constant_vars <- function(mf) {
  # don't check if columns are constant for binomial or Surv object
  mf1 <- if (NCOL(mf[, 1]) == 2 || survival::is.Surv(mf[, 1]))
    mf[, -1, drop=FALSE] else mf

  lu1 <- function(x) !all(x == 1) && length(unique(x)) == 1
  nocheck <- c("(weights)", "(offset)", "(Intercept)")
  sel <- !colnames(mf1) %in% nocheck
  is_constant <- apply(mf1[, sel, drop=FALSE], 2, lu1)
  if (any(is_constant)) {
    stop("Constant variable(s) found: ",
         paste(names(is_constant)[is_constant], collapse = ", "),
         call. = FALSE)
  }
  return(mf)
}


# Find starting
#
# @param a transition matrix
# @param a transition
match_starting <- function(t, m){
  which(apply(m, 1, function(x) return(t %in% x )))
}

# Find next state
#
# @param a transition matrix
# @param a transition
match_to <- function(t, m){
  which(apply(m, 2, function(x) return(t %in% x )) )

}

# Find competing state
#
# @param a transition matrix
# @param a initial state
# @param a ending state
match_competing <- function(m, r, t){
  m <- m[r, ][-match(names(t), names(m[r,]))]
  m[!is.na(m)]
}

# Reformulate as LHS of a formula
#
# @param x A character string or expression object
# @param as_formula Logical. If TRUE then the result is reformulated.
reformulate_lhs <- function(x) {
  #x <- deparse(x, 500L)
  x <- formula(substitute(LHS ~ 1, list(LHS = x)))
  x
}

# Extract RHS of a formula
#
# @param x A formula object
# @param as_formula Logical. If TRUE then the result is reformulated.
rhs <- function(x, as_formula = FALSE) {
  len <- length(x)
  if (len == 3L) {
    out <- x[[3L]]
  } else {
    out <- x[[2L]]
  }
  out
}

# Return a vector with valid names for elements in the list passed to the
# 'basehaz_ops' argument of a 'stan_jm' or 'stan_surv' call
#
# @param basehaz_name A character string, the type of baseline hazard.
# @return A character vector, or NA if unmatched.
get_ok_basehaz_ops <- function(basehaz_name) {
  switch(basehaz_name,
         weibull   = c(),
         bs        = c("df", "knots"),
         piecewise = c("df", "knots"),
         ms        = c("df", "knots"),
         NA)
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
handle_basehaz_surv <- function(basehaz,
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
      df <- 6L # default df for B-splines, assuming an intercept is included
    # NB this is ignored if the user specified knots

    tt <- times[status == 1] # uncensored event times
    if (is.null(knots) && !length(tt)) {
      warning2("No observed events found in the data. Censoring times will ",
               "be used to evaluate default knot locations for splines.")
      tt <- times
    }

    if (!is.null(knots)) {
      if (any(knots < min_t))
        stop2("'knots' cannot be placed before the earliest entry time.")
      if (any(knots > max_t))
        stop2("'knots' cannot be placed beyond the latest event time.")
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
      df <- 6L # default df for M-splines, assuming an intercept is included
      # NB this is ignored if the user specified knots
    }

    tt <- times[status == 1] # uncensored event times
    if (is.null(knots) && !length(tt)) {
      warning2("No observed events found in the data. Censoring times will ",
               "be used to evaluate default knot locations for splines.")
      tt <- times
    }

    if (!is.null(knots)) {
      if (any(knots < min_t))
        stop2("'knots' cannot be placed before the earliest entry time.")
      if (any(knots > max_t))
        stop2("'knots' cannot be placed beyond the latest event time.")
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

    if (!is.null(knots)) {
      if (any(knots < min_t))
        stop2("'knots' cannot be placed before the earliest entry time.")
      if (any(knots > max_t))
        stop2("'knots' cannot be placed beyond the latest event time.")
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

# Return a vector with internal knots for 'x', based on evenly spaced quantiles
#
# @param x A numeric vector.
# @param df The degrees of freedom. If specified, then 'df - degree - intercept'.
#   knots are placed at evenly spaced percentiles of 'x'. If 'iknots' is
#   specified then 'df' is ignored.
# @return A numeric vector of internal knot locations, or NULL if there are
#   no internal knots.
get_iknots <- function(x, df = 6L, degree = 3L, iknots = NULL, intercept = TRUE) {

  # obtain number of internal knots
  if (is.null(iknots)) {
    nk <- df - degree - intercept
  } else {
    nk <- length(iknots)
  }

  # validate number of internal knots
  if (nk < 0) {
    stop2("Number of internal knots cannot be negative.")
  }

  # obtain default knot locations if necessary
  if (is.null(iknots)) {
    iknots <- qtile(x, nq = nk + 1)  # evenly spaced percentiles
  }

  # return internal knot locations, ensuring they are positive
  validate_positive_scalar(iknots)

  return(iknots)
}

# Return the desired spline basis for the given knot locations
get_basis <- function(x, iknots, bknots = range(x),
                      degree = 3, intercept = TRUE,
                      type = c("bs", "is", "ms")) {
  type <- match.arg(type)
  if (type == "bs") {
    out <- splines::bs(x, knots = iknots, Boundary.knots = bknots,
                       degree = degree, intercept = intercept)
  } else if (type == "is") {
    out <- splines2::iSpline(x, knots = iknots, Boundary.knots = bknots,
                             degree = degree, intercept = intercept)
  } else if (type == "ms") {
    out <- splines2::mSpline(x, knots = iknots, Boundary.knots = bknots,
                             degree = degree, intercept = intercept)
  } else {
    stop2("'type' is not yet accommodated.")
  }
  out
}

# Return a data frame with NAs excluded
#
# @param formula The parsed model formula.
# @param data The user specified data frame.
make_model_data <- function(formula, data) {
  mf <- model.frame(formula, data, na.action = na.pass)
  include <- apply(mf, 1L, function(row) !any(is.na(row)))
  data[include, , drop = FALSE]
}

# Identify whether the type of baseline hazard requires an intercept in
# the linear predictor (NB splines incorporate the intercept into the basis).
#
# @param basehaz A list with info about the baseline hazard; see 'handle_basehaz'.
# @return A Logical.
has_intercept <- function(basehaz) {
  nm <- get_basehaz_name(basehaz)
  (nm %in% c("exp", "weibull", "gompertz"))
}

# Return the name of the baseline hazard
#
# @return A character string.
get_basehaz_name <- function(x, ind) {
  if (is.character(x))
    return(x)
  if (is.stansurv(x))
    return(x$basehaz$type_name)
  if (is.stanjm(x))
    return(x$survmod$basehaz$type_name)
  if (is.character(x$type_name))
    return(x$type_name)
  if (is.stanmstte(x))
    return(x$basehaz[[ind]]$type_name)
  stop("Bug found: could not resolve basehaz name.")
}

# Test if an object inherits a specific stanreg class
#
# @param x The object to test.
is.stanreg   <- function(x) inherits(x, "stanreg")
is.stansurv  <- function(x) inherits(x, "stansurv")
is.stanms  <- function(x) inherits(x, "stanms")
is.stanjm    <- function(x) inherits(x, "stanjm")

# Return the spline basis for the given type of baseline hazard.
#
# @param times A numeric vector of times at which to evaluate the basis.
# @param basehaz A list with info about the baseline hazard, returned by a
#   call to 'handle_basehaz'.
# @param integrate A logical, specifying whether to calculate the integral of
#   the specified basis.
# @return A matrix.
make_basis <- function(times, basehaz, integrate = FALSE) {
  N <- length(times)
  K <- basehaz$nvars
  if (!N) { # times is NULL or empty vector
    return(matrix(0, 0, K))
  }
  switch(basehaz$type_name,
         "exp"       = matrix(0, N, K), # dud matrix for Stan
         "weibull"   = matrix(0, N, K), # dud matrix for Stan
         "gompertz"  = matrix(0, N, K), # dud matrix for Stan
         "ms"        = basis_matrix(times, basis = basehaz$basis, integrate = integrate),
         "bs"        = basis_matrix(times, basis = basehaz$basis),
         "piecewise" = dummy_matrix(times, knots = basehaz$knots),
         stop2("Bug found: type of baseline hazard unknown."))
}

# Evaluate a spline basis matrix at the specified times
#
# @param time A numeric vector.
# @param basis Info on the spline basis.
# @param integrate A logical, should the integral of the basis be returned?
# @return A two-dimensional array.
basis_matrix <- function(times, basis, integrate = FALSE) {
  out <- predict(basis, times)
  if (integrate) {
    stopifnot(inherits(basis, "mSpline"))
    class(basis) <- c("matrix", "iSpline")
    out <- predict(basis, times)
  }
  aa(out)
}

# Return the fe predictor matrix for estimation
#
# @param formula The parsed model formula.
# @param model_frame The model frame.
# @return A named list with the following elements:
#   x: the fe model matrix, not centred and without intercept.
#   x_bar: the column means of the model matrix.
#   x_centered: the fe model matrix, centered.
#   N,K: number of rows (observations) and columns (predictors) in the
#     fixed effects model matrix
make_x <- function(formula, model_frame, xlevs = NULL, check_constant = TRUE) {

  # uncentred predictor matrix, without intercept
  x <- model.matrix(formula, model_frame, xlevs = xlevs)
  x <- drop_intercept(x)

  # column means of predictor matrix
  x_bar <- aa(colMeans(x))

  # centered predictor matrix
  x_centered <- sweep(x, 2, x_bar, FUN = "-")

  # identify any column of x with < 2 unique values (empty interaction levels)
  sel <- (apply(x, 2L, dplyr::n_distinct) < 2)
  if (check_constant && any(sel)) {
    cols <- paste(colnames(x)[sel], collapse = ", ")
    stop2("Cannot deal with empty interaction levels found in columns: ", cols)
  }

  nlist(x, x_centered, x_bar, N = NROW(x), K = NCOL(x))
}

# Rename the t prior as being student-t or cauchy
#
# @param prior_stuff A list with prior info returned by handle_glm_prior
# @param has A logical checking, for example, whether the model has_predictors,
#   has_intercept, has_assoc, etc
rename_t_and_cauchy <- function(prior_stuff, has) {
  if (has && prior_stuff$prior_dist_name %in% "t") {
    if (all(prior_stuff$prior_df == 1)) {
      prior_stuff$prior_dist_name <- "cauchy"
    } else {
      prior_stuff$prior_dist_name <- "student_t"
    }
  }
  return(prior_stuff)
}

# Get name of auxiliary parameters for event submodel
#
# @param basehaz A list with information about the baseline hazard
.rename_e_aux <- function(basehaz) {
  nm <- basehaz$type_name
  switch(nm,
         weibull   = "weibull-shape",
         gompertz  = "gompertz-scale",
         bs        = "B-spline-coefficients",
         ms        = "M-spline-coefficients",
         piecewise = "piecewise-coefficients",
         NA)
}



# Set arguments for sampling
#
# Prepare a list of arguments to use with \code{rstan::sampling} via
# \code{do.call}.
#
# @param object The stanfit object to use for sampling.
# @param user_dots The contents of \code{...} from the user's call to
#   the \code{stan_*} modeling function.
# @param user_adapt_delta The value for \code{adapt_delta} specified by the
#   user.
# @param prior Prior distribution list (can be NULL).
# @param ... Other arguments to \code{\link[rstan]{sampling}} not coming from
#   \code{user_dots} (e.g. \code{data}, \code{pars}, \code{init}, etc.)
# @return A list of arguments to use for the \code{args} argument for
#   \code{do.call(sampling, args)}.
set_sampling_args <- function(object, prior, user_dots = list(),
                              user_adapt_delta = NULL, ...) {
  args <- list(object = object, ...)
  unms <- names(user_dots)
  for (j in seq_along(user_dots)) {
    args[[unms[j]]] <- user_dots[[j]]
  }
  defaults <- default_stan_control(prior = prior,
                                   adapt_delta = user_adapt_delta)

  if (!"control" %in% unms) {
    # no user-specified 'control' argument
    args$control <- defaults
  } else {
    # user specifies a 'control' argument
    if (!is.null(user_adapt_delta)) {
      # if user specified adapt_delta argument to stan_* then
      # set control$adapt_delta to user-specified value
      args$control$adapt_delta <- user_adapt_delta
    } else {
      # use default adapt_delta for the user's chosen prior
      args$control$adapt_delta <- defaults$adapt_delta
    }
    if (is.null(args$control$max_treedepth)) {
      # if user's 'control' has no max_treedepth set it to rstanarm default
      args$control$max_treedepth <- defaults$max_treedepth
    }
  }
  args$save_warmup <- FALSE

  return(args)
}

# Default control arguments for sampling
#
# Called by set_sampling_args to set the default 'control' argument for
# \code{rstan::sampling} if none specified by user. This allows the value of
# \code{adapt_delta} to depend on the prior.
#
# @param prior Prior distribution list (can be NULL).
# @param adapt_delta User's \code{adapt_delta} argument.
# @param max_treedepth Default for \code{max_treedepth}.
# @return A list with \code{adapt_delta} and \code{max_treedepth}.
default_stan_control <- function(prior, adapt_delta = NULL,
                                 max_treedepth = 15L) {
  if (!length(prior)) {
    if (is.null(adapt_delta)) adapt_delta <- 0.95
  } else if (is.null(adapt_delta)) {
    adapt_delta <- switch(prior$dist,
                          "R2" = 0.99,
                          "hs" = 0.99,
                          "hs_plus" = 0.99,
                          "lasso" = 0.99,
                          "product_normal" = 0.99,
                          0.95) # default
  }
  nlist(adapt_delta, max_treedepth)
}
