#------------------------------
# Below are code chunks taken from the 'rstanarm' R package, obtained
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

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

# Unlist the result from an lapply call with recursive = FALSE
#
# @param X,FUN,... Same as lapply
nonrecuapply <- function(X, FUN, ...) {
  unlist(lapply(X, FUN, ...), recursive = FALSE)
}

# A refactored version of mapply with SIMPLIFY = FALSE
#
# @param FUN,... Same as mapply
# @param arg Passed to MoreArgs
xapply <- function(..., FUN, args = NULL) {
  mapply(FUN, ..., MoreArgs = args, SIMPLIFY = FALSE)
}


# Stop without printing call
stop2    <- function(...) stop(..., call. = FALSE)

# Immediate warning without printing call
warning2 <- function(...) warning(..., immediate. = TRUE, call. = FALSE)

# Shorthand for suppress warnings
SW <- function(expr) base::suppressWarnings(expr)
SM <- function(expr) base::suppressMessages(exp)
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
uu <- function(x, ...) unique(unlist(x, ...))
# Safe deparse
safe_deparse <- function(expr) deparse(expr, 500L)

# Transpose function that can handle NULL objects
#
# @param x A matrix, a vector, or otherwise (e.g. NULL)
transpose <- function(x) {
  if (is.matrix(x) || is.vector(x)) {
    t(x)
  } else {
    x
  }
}


# If x is NULL then return an empty object of the specified 'type'
#
# @param x An object to test whether it is null.
# @param type The type of empty object to return if x is null.
convert_null <- function(x, type = c("double", "integer", "matrix",
                                     "arraydouble", "arrayinteger")) {
  if (!is.null(x)) {
    return(x)
  } else if (type == "double") {
    return(double(0))
  } else if (type == "integer") {
    return(integer(0))
  } else if (type == "matrix") {
    return(matrix(0,0,0))
  } else if (type == "arraydouble") {
    return(as.array(double(0)))
  } else if (type == "arrayinteger") {
    return(as.array(integer(0)))
  } else {
    stop("Input type not valid.")
  }
}
# Translate group/factor IDs into integer values
#
# @param x A vector of group/factor IDs
groups <- function(x) {
  if (!is.null(x)) {
    as.integer(as.factor(x))
  } else {
    x
  }
}

# If a is NULL (and Inf, respectively) return b, otherwise just return a
# @param a,b Objects
`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}
`%ORifINF%` <- function(a, b) {
  if (a == Inf) b else a
}
`%!in%` <- function(x,y)!('%in%'(x,y))

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

# Return the names for the coefficients
get_beta_name_ymod <- function(x) {
  nms <- colnames(x$x$xtemp)
  if (!is.null(nms)) paste0(x$stub, "|", nms) else NULL
}

get_aux_name_ymod <- function(x, ...) {
  switch(x$family$family,
         gaussian         = paste0(x$stub, "|sigma"),
         Gamma            = paste0(x$stub, "|shape"),
         inverse.gaussian = paste0(x$stub, "|lambda"),
         neg_binomial_2   = paste0(x$stub, "|reciprocal_dispersion"),
         NULL)
}

get_int_name_ymod <- function(x, ...) {
  if (x$intercept_type$number) paste0(x$stub, "|(Intercept)") else NULL
}

# Return the name for the mean_PPD
get_ppd_name <- function(x, ...) {
  paste0(x$stub, "|mean_PPD")
}

# Add the variables in ...'s to the RHS of a model formula
#
# @param x A model formula.
# @param ... Character strings, the variable names.
addto_formula <- function(x, ...) {
  rhs_terms   <- terms(reformulate_rhs(rhs(x)))
  intercept   <- attr(rhs_terms, "intercept")
  term_labels <- attr(rhs_terms, "term.labels")
  reformulate(c(term_labels, c(...)), response = lhs(x), intercept = intercept)
}

# Get the posterior sample size
#
# @param x A stanreg object
# @return the posterior sample size (or size of sample from approximate posterior)
posterior_sample_size <- function(x) {
  validate_stanmstte_object(x)
  if (used.optimizing(x)) {
    return(NROW(x$asymptotic_sampling_dist))
  }
  pss <- x$stanfit@sim$n_save
  if (used.variational(x))
    return(pss)
  sum(pss - x$stanfit@sim$warmup2)
}

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

# Consistent error message when binomial models with greater than
# one trial are not allowed
#
STOP_binomial <- function() {
  stop2("Binomial models with number of trials greater than one ",
        "are not allowed (i.e. only bernoulli models are allowed).")
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

# Verify that outcome values match support implied by family object
#
# @param y outcome variable
# @param family family object
# @return y (possibly slightly modified) unless an error is thrown
#
validate_glm_outcome_support <- function(y, family) {
  .is_count <- function(x) {
    all(x >= 0) && all(abs(x - round(x)) < .Machine$double.eps^0.5)
  }

  fam <- family$family

  if (!is.binomial(fam)) {
    # make sure y has ok dimensions (matrix only allowed for binomial models)
    if (length(dim(y)) > 1) {
      if (NCOL(y) == 1) {
        y <- y[, 1]
      } else {
        stop("Except for binomial models the outcome variable ",
             "should not have multiple columns.",
             call. = FALSE)
      }
    }

    # check that values match support for non-binomial models
    if (is.gaussian(fam)) {
      return(y)
    } else if (is.gamma(fam) && any(y <= 0)) {
      stop("All outcome values must be positive for gamma models.",
           call. = FALSE)
    } else if (is.ig(fam) && any(y <= 0)) {
      stop("All outcome values must be positive for inverse-Gaussian models.",
           call. = FALSE)
    } else if (is.poisson(fam) && !.is_count(y)) {
      stop("All outcome values must be counts for Poisson models",
           call. = FALSE)
    } else if (is.nb(fam) && !.is_count(y)) {
      stop("All outcome values must be counts for negative binomial models",
           call. = FALSE)
    }
  } else { # binomial models
    if (NCOL(y) == 1L) {
      if (is.numeric(y) || is.logical(y))
        y <- as.integer(y)
      if (is.factor(y))
        y <- fac2bin(y)
      if (!all(y %in% c(0L, 1L)))
        stop("All outcome values must be 0 or 1 for Bernoulli models.",
             call. = FALSE)
    } else if (isTRUE(NCOL(y) == 2L)) {
      if (!.is_count(y))
        stop("All outcome values must be counts for binomial models.",
             call. = FALSE)
    } else {
      stop("For binomial models the outcome should be a vector or ",
           "a matrix with 2 columns.",
           call. = FALSE)
    }
  }

  return(y)
}

# Check family argument
#
# @param f The \code{family} argument specified by user (or the default).
# @return If no error is thrown, then either \code{f} itself is returned (if
#   already a family) or the family object created from \code{f} is returned (if
#   \code{f} is a string or function).
validate_family <- function(f) {
  if (is.character(f))
    f <- get(f, mode = "function", envir = parent.frame(2))
  if (is.function(f))
    f <- f()
  if (!is(f, "family"))
    stop("'family' must be a family.", call. = FALSE)

  return(f)
}

# Check if x and any objects in ... were all NULL or not
#
# @param x The first object to use in the comparison
# @param ... Any additional objects to include in the comparison
# @param error If TRUE then return an error if all objects aren't
#   equal with regard to the 'is.null' test.
# @return If error = TRUE, then an error if all objects aren't
#   equal with regard to the 'is.null' test. Otherwise, a logical
#   specifying whether all objects were equal with regard to the
#   'is.null' test.
supplied_together <- function(x, ..., error = FALSE) {
  dots <- list(...)
  for (i in dots) {
    if (!identical(is.null(x), is.null(i))) {
      if (error) {
        nm_x <- deparse(substitute(x))
        nm_i <- deparse(substitute(i))
        stop2(nm_x, " and ", nm_i, " must be supplied together.")
      } else {
        return(FALSE) # not supplied together, ie. one NULL and one not NULL
      }
    }
  }
  return(TRUE) # supplied together, ie. all NULL or all not NULL
}


# Throw error if any transition is not present in the dataset.
#
# @param x The long-format dataframe to test.
# @param n number of transition
validate_n_trans <- function(x, n){
  if(!all(unique(x$trans) %in% seq_len(n)))
    stop2("Attempted to estimate more transitions than found in the dataset, reconsider the model given your dataset.")
}

# Function to check if the response vector is real or integer
#
# @param family A family object
# @return A logical specify whether the response is real (TRUE) or integer (FALSE)
check_response_real <- function(family) {
  !(family$family %in% c("binomial", "poisson", "neg_binomial_2"))
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

# Return a list (or vector if unlist = TRUE) which
# contains the embedded elements in list x named y
fetch <- function(x, y, z = NULL, zz = NULL, null_to_zero = FALSE,
                  pad_length = NULL, unlist = FALSE) {
  ret <- lapply(x, `[[`, y)
  if (!is.null(z))
    ret <- lapply(ret, `[[`, z)
  if (!is.null(zz))
    ret <- lapply(ret, `[[`, zz)
  if (null_to_zero)
    ret <- lapply(ret, function(i) ifelse(is.null(i), 0L, i))
  if (!is.null(pad_length)) {
    padding <- rep(list(0L), pad_length - length(ret))
    ret <- c(ret, padding)
  }
  if (unlist) unlist(ret) else ret
}

# Wrapper for using fetch with unlist = TRUE and
# returning array. Also converts logical to integer.
fetch_array <- function(x, y, z = NULL, zz = NULL, null_to_zero = FALSE,
                        pad_length = NULL) {
  val <- fetch(x = x, y = y, z = z, zz = zz, null_to_zero = null_to_zero,
               pad_length = pad_length, unlist = TRUE)
  if (is.logical(val))
    val <- as.integer(val)
  as.array(val)
}

# Wrapper for using fetch with unlist = TRUE
fetch_ <- function(x, y, z = NULL, zz = NULL, null_to_zero = FALSE,
                   pad_length = NULL) {
  fetch(x = x, y = y, z = z, zz = zz, null_to_zero = null_to_zero,
        pad_length = pad_length, unlist = TRUE)
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

# Test for a given family
#
# @param x A character vector (probably x = family(fit)$family)
is.binomial <- function(x) x == "binomial"
is.gaussian <- function(x) x == "gaussian"
is.gamma    <- function(x) x == "Gamma"
is.ig       <- function(x) x == "inverse.gaussian"
is.nb       <- function(x) x == "neg_binomial_2"
is.poisson  <- function(x) x == "poisson"
is.beta     <- function(x) x == "beta" || x == "Beta regression"

# Test for a given estimation method
#
# @param x A stanreg object.
used.optimizing <- function(x) {
  x$algorithm == "optimizing"
}
used.sampling <- function(x) {
  x$algorithm == "sampling"
}
used.variational <- function(x) {
  x$algorithm %in% c("meanfield", "fullrank")
}

# Invert 'is.null'
not.null <- function(x) { !is.null(x) }

# loo internal ----------------------------------------------------------------
validate_k_threshold <- function(k) {
  if (!is.numeric(k) || length(k) != 1) {
    stop("'k_threshold' must be a single numeric value.",
         call. = FALSE)
  } else if (k < 0) {
    stop("'k_threshold' < 0 not allowed.",
         call. = FALSE)
  } else if (k > 1) {
    warning(
      "Setting 'k_threshold' > 1 is not recommended.",
      "\nFor details see the PSIS-LOO section in help('loo-package', 'loo').",
      call. = FALSE
    )
  }
}


# Split a vector or matrix into a specified number of segments and return
# each segment as an element of a list. The matrix method allows splitting
# across the column (bycol = TRUE) or row margin (bycol = FALSE).
#
# @param x A vector or matrix.
# @param n_segments Integer specifying the number of segments.
# @param bycol Logical, should a matrix be split along the column or row margin?
# @return A list with n_segments elements.
split2 <- function(x, n_segments = 1, ...) {
  UseMethod("split2")
}

split2.vector <- function(x, n_segments = 1, ...) {
  len <- length(x)
  segment_length <- len %/% n_segments
  if (!len == (segment_length * n_segments))
    stop("Dividing x by n_segments does not result in an integer.")
  split(x, rep(1:n_segments, each = segment_length))
}

split2.matrix <- function(x, n_segments = 1, bycol = TRUE) {
  len <- if (bycol) ncol(x) else nrow(x)
  segment_length <- len %/% n_segments
  if (!len == (segment_length * n_segments))
    stop("Dividing x by n_segments does not result in an integer.")
  lapply(1:n_segments, function(k) {
    if (bycol) x[, (k-1) * segment_length + 1:segment_length, drop = FALSE] else
      x[(k-1) * segment_length + 1:segment_length, , drop = FALSE]})
}

# Split a vector or matrix into a specified number of segments
# (see rstanarm:::split2) and then reduce it using 'FUN'
split_and_reduce <- function(x, n_segments = 1, bycol = TRUE, FUN = '+') {
  splitted_x <- split2(x, n_segments = n_segments, bycol = bycol)
  Reduce(FUN, splitted_x)
}

# chain_id to pass to loo::relative_eff
chain_id_for_loo <- function(object) {
  dims <- dim(object$stanfit)[1:2]
  n_iter <- dims[1]
  n_chain <- dims[2]
  rep(1:n_chain, each = n_iter)
}


# check if discrete or continuous
# @param object stanreg object
is_discrete <- function(object) {
  if (inherits(object, "polr"))
    return(TRUE)
  if (inherits(object, "stansurv"))
    return(FALSE)
  if (inherits(object, "stanmstte"))
    return(FALSE)
  if (inherits(object, "stanmvreg")) {
    fams <- fetch(family(object), "family")
    res <- sapply(fams, function(x)
      is.binomial(x) || is.poisson(x) || is.nb(x))
    return(res)
  }
  fam <- family(object)$family
  is.binomial(fam) || is.poisson(fam) || is.nb(fam)
}

# Calculate a SHA1 hash of y
# @param x stanreg object
# @param ... Passed to digest::sha1
#
hash_y <- function(x, ...) {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("Please install the 'digest' package.")
  validate_stanmstte_object(x)
  y <- get_y(x)
  for(i in seq_along(y)){
    attributes( y[[i]] ) <- NULL
  }
  digest::sha1(x = y, ...)
}

# validate objects for model comparison
validate_loos <- function(loos = list()) {
  if (length(loos) <= 1)
    stop("At least two objects are required for model comparison.",
         call. = FALSE)

  is_loo <- sapply(loos, is.loo)
  is_waic <- sapply(loos, is.waic)
  is_kfold <- sapply(loos, is.kfold)
  if (!all(is_loo))
    stop("All objects must have class 'loo'", call. = FALSE)
  if ((any(is_waic) && !all(is_waic) ||
       (any(is_kfold) && !all(is_kfold))))
    stop("Can't mix objects computed using 'loo', 'waic', and 'kfold'.",
         call. = FALSE)

  yhash <- lapply(loos, attr, which = "yhash")
  yhash_check <- sapply(yhash, function(x) {
    isTRUE(all.equal(x, yhash[[1]]))
  })
  if (!all(yhash_check))
    stop("Not all models have the same y variable.", call. = FALSE)

  discrete <- sapply(loos, attr, which = "discrete")
  if (!all(discrete == discrete[1]))
    stop("Discrete and continuous observation models can't be compared.",
         call. = FALSE)

  setNames(loos, nm = lapply(loos, attr, which = "name"))
}

# model formula to store in loo object
# @param x stanreg object
loo_model_formula <- function(x) {
  form <- try(formula(x), silent = TRUE)
  if (inherits(form, "try-error") || is.null(form)) {
    form <- "formula not found"
  }
  return(form)
}


# ------------- Helpers ---------------#
# Count the number of unique values
#
# @param x A vector or list
n_distinct <- function(x) {
  length(unique(x))
}
# Left join NA to 0 for sparse matrix
full_join_NA <- function(x, y, ...) {
  dplyr::full_join(x = x, y = y, by = ...) %>%
    dplyr::mutate_each(dplyr::funs(replace(., which(is.na(.)), 0)))
}
# Combine pars and regex_pars
#
# @param x stanreg object
# @param pars Character vector of parameter names
# @param regex_pars Character vector of patterns
collect_pars <- function(x, pars = NULL, regex_pars = NULL) {
  if (is.null(pars) && is.null(regex_pars))
    return(NULL)
  if (!is.null(pars))
    pars[pars == "varying"] <- "b"
  if (!is.null(regex_pars))
    pars <- c(pars, grep_for_pars(x, regex_pars))
  unique(pars)
}
# Grep for "b" parameters (ranef)
#
# @param x Character vector (often rownames(fit$stan_summary))
# @param ... Passed to grep
b_names <- function(x, ...) {
  grep("^b\\[", x, ...)
}
# Extract parameters from stanmat and return as a list
#
# @param object A stanmvreg or stansurv object
# @param stanmat A matrix of posterior draws, may be provided if the desired
#   stanmat is only a subset of the draws from as.matrix(object$stanfit)
# @return A named list
extract_pars <- function(object, ...) {
  UseMethod("extract_pars")
}

extract_pars.stanmstte <- function(object, stanmat = NULL, means = FALSE) {
  validate_stanmstte_object(object)
  if (is.null(stanmat))
    stanmat <- as.matrix(object$stanfit)
  if (means)
    stanmat <- t(colMeans(stanmat)) # return posterior means
  nms_beta <- lapply(object$x, function(x) colnames(x) )
  nms_beta <- uapply(seq_len(object$n_trans), function(i)
     append_trans_to_pars(nms_beta[[i]], i, object$transition_labels[i])
  )
  nms_tde  <- get_smooth_name(object$s_cpts, type = "smooth_coefs")
  nms_smth <- get_smooth_name(object$s_cpts, type = "smooth_sd")
  nms_int  <- lapply(object$basehaz, function(b) get_int_name_basehaz(b) )
  nms_int <- uapply(seq_len(object$n_trans), function(i)
    append_trans_to_pars(nms_int[[i]], i, object$transition_labels[i]) )
  nms_aux  <- lapply(object$basehaz, function(b)  get_aux_name_basehaz(b) )
  nms_aux <- uapply(seq_len(object$n_trans), function(i)
    append_trans_to_pars(nms_aux[[i]], i, object$transition_labels[i]) )
  alpha    <- stanmat[, nms_int,  drop = FALSE]
  beta     <- stanmat[, nms_beta, drop = FALSE]
  beta_tde <- stanmat[, nms_tde,  drop = FALSE]
  aux      <- stanmat[, nms_aux,  drop = FALSE]
  smooth   <- stanmat[, nms_smth, drop = FALSE]
  nlist(alpha, beta, beta_tde, aux, smooth, stanmat)
}

# Return the name of the tde spline coefs or smoothing parameters.
#
# @param x The predictor matrix for the time-dependent effects, with column names.
# @param type The type of information about the smoothing parameters to return.
# @return A character or numeric vector, depending on 'type'.
get_smooth_name <- function(x, type = "smooth_coefs") {

  if (is.null(x) || !ncol(x))
    return(NULL)

  nms <- gsub(":bs\\(times__.*\\)[0-9]*$", "", colnames(x))
  tally   <- table(nms)
  indices <- uapply(tally, seq_len)
  suffix  <- paste0(":tde-spline-coef", indices)

  switch(type,
         smooth_coefs = paste0(nms, suffix),
         smooth_sd    = paste0("smooth_sd[", unique(nms), "]"),
         smooth_map   = rep(seq_along(tally), tally),
         smooth_vars  = unique(nms),
         stop2("Bug found: invalid input to 'type' argument."))
}

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

append_trans_to_pars <- function(x, i, labs){
  if(is.null(x) ) return(NULL)
  if(is.null(labs)){
    paste0(x," trans(",i,")" )
  } else {
    paste0(x," trans(",labs,")" )
  }
}

# Append a string (prefix) to the column names of a matrix or array
append_prefix_to_colnames <- function(x, str) {
  if (ncol(x)) set_colnames(x, paste0(str, colnames(x))) else x
}

set_colnames <- function(x, names) { colnames(x) <- names; x }

append_title <- function(object){
  if(is.null(object$transition_labels)){
    paste("\nTransition", seq_len(object$n_trans),"\n")
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
  if(!missing(ind))
    object[["y"]][[ind]] %ORifNULL% model.response(model.frame(object))
  else {
      lapply(object$model_frame, function(m) model.response(m) )
  }
}
#' @export
get_x.default <- function(object, ind, ...) {
  object[["x"]][[ind]] %ORifNULL% model.matrix(object)
}

# Methods for creating linear predictor
#
# Make linear predictor vector from x and point estimates for beta, or linear
# predictor matrix from x and full posterior sample of beta.
#
# @param beta A vector or matrix or parameter estimates.
# @param x Predictor matrix.
# @param offset Optional offset vector.
# @return A vector or matrix.
linear_predictor <- function(beta, x, offset = NULL) {
  UseMethod("linear_predictor")
}
linear_predictor.default <- function(beta, x, offset = NULL) {
  eta <- as.vector(if (NCOL(x) == 1L) x * beta else x %*% beta)
  if (length(offset))
    eta <- eta + offset

  return(eta)
}
linear_predictor.matrix <- function(beta, x, offset = NULL) {
  if (NCOL(beta) == 1L)
    beta <- as.matrix(beta)
  eta <- beta %*% t(x)
  if (length(offset))
    eta <- sweep(eta, 2L, offset, `+`)

  return(eta)
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

# Throw error if object isn't a stanreg object
#
# @param x The object to test.
validate_stanmstte_object <- function(x, call. = FALSE) {
  if (!is.stanmstte(x))
    stop("Object is not a stanmstee object.", call. = call.)
}

# Validate newdataLong and newdataEvent arguments
#
# @param object A stanmvreg object
# @param newdataLong A data frame, or a list of data frames
# @param newdataEvent A data frame
# @param duplicate_ok A logical. If FALSE then only one row per individual is
#   allowed in the newdataEvent data frame
# @param response A logical specifying whether the longitudinal response
#   variable must be included in the new data frame
# @return A list of validated data frames
validate_newdatas <- function(object, newdataLong = NULL, newdataEvent = NULL,
                              duplicate_ok = FALSE, response = TRUE) {
  validate_stanmstte_object(object)
  id_var <- object$id_var
  newdatas <- list()
  if (!is.null(newdataLong)) {
    if (!is(newdataLong, "list"))
      newdataLong <- rep(list(newdataLong), get_M(object))
    dfcheck <- sapply(newdataLong, is.data.frame)
    if (!all(dfcheck))
      stop("'newdataLong' must be a data frame or list of data frames.", call. = FALSE)
    nacheck <- sapply(seq_along(newdataLong), function(m) {
      if (response) { # newdataLong needs the reponse variable
        fmL <- formula(object, m = m)
      } else { # newdataLong only needs the covariates
        fmL <- formula(object, m = m)[c(1,3)]
      }
      all(!is.na(get_all_vars(fmL, newdataLong[[m]])))
    })
    if (!all(nacheck))
      stop("'newdataLong' cannot contain NAs.", call. = FALSE)
    newdatas <- c(newdatas, newdataLong)
  }
  if (!is.null(newdataEvent)) {
    if (!is.data.frame(newdataEvent))
      stop("'newdataEvent' must be a data frame.", call. = FALSE)
    if (response) { # newdataEvent needs the reponse variable
      fmE <- formula(object, m = "Event")
    } else { # newdataEvent only needs the covariates
      fmE <- formula(object, m = "Event")[c(1,3)]
    }
    dat <- get_all_vars(fmE, newdataEvent)
    dat[[id_var]] <- newdataEvent[[id_var]] # include ID variable in event data
    if (any(is.na(dat)))
      stop("'newdataEvent' cannot contain NAs.", call. = FALSE)
    if (!duplicate_ok && any(duplicated(newdataEvent[[id_var]])))
      stop("'newdataEvent' should only contain one row per individual, since ",
           "time varying covariates are not allowed in the prediction data.")
    newdatas <- c(newdatas, list(Event = newdataEvent))
  }
  if (length(newdatas)) {
    idvar_check <- sapply(newdatas, function(x) id_var %in% colnames(x))
    if (!all(idvar_check))
      STOP_no_var(id_var)
    ids <- lapply(newdatas, function(x) unique(x[[id_var]]))
    sorted_ids <- lapply(ids, sort)
    if (!length(unique(sorted_ids)) == 1L)
      stop("The same subject ids should appear in each new data frame.")
    if (!length(unique(ids)) == 1L)
      stop("The subject ids should be ordered the same in each new data frame.")
    return(newdatas)
  } else return(NULL)
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

# Return the number of longitudinal submodels
#
# @param object A stanmvreg object
get_T <- function(object) {
  validate_stanmvreg_object(object)
  return(object$n_trans)
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

get_id_var <- function(x){
  if(!is.null(x$id_var)){
    x$id_var
  } else {
    as.character( rownames(x) )
  }
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

# Check if a fitted model (stanreg object) has weights
#
# @param x stanreg object
# @return Logical. Only TRUE if x$weights has positive length and the elements
#   of x$weights are not all the same.
#
model_has_weights <- function(x) {
  wts <- x[["weights"]]
  if (!length(wts)) {
    FALSE
  } else if (all(wts == wts[1])) {
    FALSE
  } else {
    TRUE
  }
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
has_intercept <- function(basehaz, ind = NULL) {
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

# Convert a standardised quadrature node to an unstandardised value based on
# the specified integral limits
#
# @param x An unstandardised quadrature node
# @param a The lower limit(s) of the integral, possibly a vector
# @param b The upper limit(s) of the integral, possibly a vector
unstandardise_qpts <- function(x, a, b) {
  if (!identical(length(x), 1L) || !is.numeric(x))
    stop("'x' should be a single numeric value.", call. = FALSE)
  if (!all(is.numeric(a), is.numeric(b)))
    stop("'a' and 'b' should be numeric.", call. = FALSE)
  if (!length(a) %in% c(1L, length(b)))
    stop("'a' and 'b' should be vectors of length 1, or, be the same length.", call. = FALSE)
  if (any((b - a) < 0))
    stop("The upper limits for the integral ('b' values) should be greater than ",
         "the corresponding lower limits for the integral ('a' values).", call. = FALSE)
  ((b - a) / 2) * x + ((b + a) / 2)
}

# Convert a standardised quadrature weight to an unstandardised value based on
# the specified integral limits
#
# @param x An unstandardised quadrature weight
# @param a The lower limit(s) of the integral, possibly a vector
# @param b The upper limit(s) of the integral, possibly a vector
unstandardise_qwts <- function(x, a, b) {
  if (!identical(length(x), 1L) || !is.numeric(x))
    stop("'x' should be a single numeric value.", call. = FALSE)
  if (!all(is.numeric(a), is.numeric(b)))
    stop("'a' and 'b' should be numeric.", call. = FALSE)
  if (!length(a) %in% c(1L, length(b)))
    stop("'a' and 'b' should be vectors of length 1, or, be the same length.", call. = FALSE)
  if (any((b - a) < 0))
    stop("The upper limits for the integral ('b' values) should be greater than ",
         "the corresponding lower limits for the integral ('a' values).", call. = FALSE)
  ((b - a) / 2) * x
}


# Replicate rows of a matrix or data frame
#
# @param x A matrix or data frame.
# @param ... Arguments passed to 'rep', namely 'each' or 'times'.
rep_rows <- function(x, ...) {
  if (is.null(x) || !nrow(x)) {
    return(x)
  } else if (is.matrix(x) || is.data.frame(x)) {
    x <- x[rep(1:nrow(x), ...), , drop = FALSE]
  } else {
    stop2("'x' must be a matrix or data frame.")
  }
  x
}


# Promote a character variable to a factor
#
# @param x The variable to potentially promote
promote_to_factor <- function(x) {
  if (is.character(x)) as.factor(x) else x
}
