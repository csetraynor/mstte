#' Print method for stanmstte objects
#'
#' The \code{print} method for stanmstte objects displays a compact summary of the
#' fitted model. See the \strong{Details} section below for descriptions of the
#' different components of the printed output. For additional summary statistics
#' and diagnostics use the \code{\link[=summary.stanmstte]{summary}} method.
#'
#' @export
#' @method print stanmstte
#' @param digits Number of digits to use for formatting numbers.
#' @param ... Ignored.
#' @return Returns \code{x}, invisibly.
#' @details
#' \subsection{Point estimates}{
#' Regardless of the estimation algorithm, point estimates are medians computed
#' from simulations. For models fit using MCMC (\code{"sampling"}) the posterior
#' sample is used. For the \code{"meanfield"} and \code{"fullrank"} variational
#' approximations, draws from the variational approximation to the posterior are
#' used. In all cases, the point estimates reported are the same as the values
#' returned by \code{\link[=coef.stanmstte]{coef}}.
#' }
#' \subsection{Uncertainty estimates (MAD_SD)}{
#' The standard deviations reported (labeled \code{MAD_SD} in the print output)
#' are computed from the same set of draws described above and are proportional
#' to the median absolute deviation (\code{\link[stats]{mad}}) from the median.
#' Compared to the raw posterior standard deviation, the MAD_SD will be
#' more robust for long-tailed distributions. These are the same as the values
#' returned by \code{\link[=se.stanmstte]{se}}.
#' }
#' \subsection{Additional output}{
#' \itemize{
#' \item The median and MAD_SD are also reported for \code{mean_PPD}, the sample
#' average posterior predictive distribution of the outcome. This is useful as a
#' quick diagnostic. A useful heuristic is to check if \code{mean_PPD} is
#' plausible when compared to \code{mean(y)}. If it is plausible then this does
#' \item For multi-state model the estimates are presented seprated for each transtion hazard
#' \item For joint longitudinal and time-to-event (see \code{\link{stan_jm}}) models
#' the estimates are presented separately for each of the distinct submodels.
#' }
#' }
#'
#'
print.stanmstte <- function(x, digits = 3, ...) {
  ##
  cat("\nAnalysis of independent multi-state risks data \n")
  title <- append_title(x)

  for(i in seq_along(x$basehaz)){
    cat("\n ---------------- \n")
    cat(title[[i]])
    cat("\n baseline hazard:", basehaz_string(x$basehaz[[i]]) )
    cat("\n formula:        ", formula_string(formula(x)[[i]]$formula) )
    cat("\n observations:   ", x$nobs[[i]])
    cat("\n events:         ", x$nevents[[i]], percent_string(x$nevents[[i]], x$nobs[[i]]))
    if (x$nlcens[[i]] > 0)
      cat("\n left censored:  ", x$nlcens[[i]], percent_string(x$nlcens[[i]], x$nobs[[i]]))
    if (x$nrcens[[i]] > 0)
      cat("\n right censored: ", x$nrcens[[i]], percent_string(x$nrcens[[i]], x$nobs[[i]]))
    if (x$nicens[[i]] > 0)
      cat("\n interval cens.: ", x$nicens[[i]], percent_string(x$nicens[[i]], x$nobs[[i]]))
    cat("\n delayed entry:  ", yes_no_string(x$ndelayed[[i]]))

    cat("\nEstimates:\n")
    mat <- as.matrix(x$stanfit) # don't used as.matrix.stanreg method b/c want access to mean_PPD
    aux_nms <- .aux_name(x)
    nms <- setdiff(rownames(x$stan_summary), c("log-posterior", aux_nms))
    nms <- grep(get_transition_name(x, i), nms, value = TRUE, fixed = TRUE)
    ppd_nms <- grep("^mean_PPD", nms, value = TRUE)
    nms <- setdiff(nms, ppd_nms)
    coef_mat <- mat[, nms, drop = FALSE]
    estimates <- .median_and_madsd(coef_mat)
    if(has_intercept(x$basehaz[[i]])){
      nms_int  <- append_trans( get_int_name_basehaz(get_basehaz(x, i)), i, x$transition_labels[i])
    } else {
      nms_int = NULL
    }

    if(get_basehaz_name(x$basehaz[[i]]) != "exp"){
      nms_aux  <- append_trans(get_aux_name_basehaz(get_basehaz(x, i)), i, x$transition_labels[i])
    } else {
      nms_aux = NULL
    }
    if(ncol(get_x(x, i)) > 0){
      nms_beta <-  append_trans(colnames(get_x(x, i)), i, x$transition_labels[i])
    } else {
      nms_beta <- NULL
    }

    estimates <- cbind(estimates,
            "exp(Median)" = c(rep(NA, length(nms_int)), exp(estimates[nms_beta, "Median"]), rep(NA, length(nms_aux)))
      )
    .printfr(estimates, digits, ...)

  }

  cat("\n ---------------- \n")
  cat("* For help interpreting the printed output see ?print.stanidm\n")
  cat("* For info on the priors used see ?prior_summary.stanreg\n")
}

#' Summary method for stanmstte objects
#'
#' Summaries of parameter estimates and MCMC convergence diagnostics
#' (Monte Carlo error, effective sample size, Rhat).
#'
#' @export
#' @method summary stanmstte
#'
#' @param ... Currently ignored.
#' @param pars An optional character vector specifying a subset of parameters to
#'   display. Parameters can be specified by name or several shortcuts can be
#'   used. Using \code{pars="beta"} will restrict the displayed parameters to
#'   only the regression coefficients (without the intercept). \code{"alpha"}
#'   can also be used as a shortcut for \code{"(Intercept)"}. If the model has
#'   varying intercepts and/or slopes they can be selected using \code{pars =
#'   "varying"}.
#'
#'   In addition, for \code{stanmvreg} objects there are some additional shortcuts
#'   available. Using \code{pars = "long"} will display the
#'   parameter estimates for the longitudinal submodels only (excluding group-specific
#'   pparameters, but including auxiliary parameters).
#'   Using \code{pars = "event"} will display the
#'   parameter estimates for the event submodel only, including any association
#'   parameters.
#'   Using \code{pars = "assoc"} will display only the
#'   association parameters.
#'   Using \code{pars = "fixef"} will display all fixed effects, but not
#'   the random effects or the auxiliary parameters.
#'    \code{pars} and \code{regex_pars} are set to \code{NULL} then all
#'   fixed effect regression coefficients are selected, as well as any
#'   auxiliary parameters and the log posterior.
#'
#'   If \code{pars} is \code{NULL} all parameters are selected for a \code{stanreg}
#'   object, while for a \code{stanmvreg} object all
#'   fixed effect regression coefficients are selected as well as any
#'   auxiliary parameters and the log posterior. See
#'   \strong{Examples}.
#' @param probs For models fit using MCMC or one of the variational algorithms,
#'   an optional numeric vector of probabilities passed to
#'   \code{\link[stats]{quantile}}.
#' @param digits Number of digits to use for formatting numbers when printing.
#'   When calling \code{summary}, the value of digits is stored as the
#'   \code{"print.digits"} attribute of the returned object.
#'
#' @return The \code{summary} method returns an object of class
#'   \code{"summary.stanmstte"} (or \code{"summary.stanmvreg"}, inheriting
#'   \code{"summary.stanmstte"}), which is a matrix of
#'   summary statistics and
#'   diagnostics, with attributes storing information for use by the
#'   \code{print} method. The \code{print} method for \code{summary.stanmstte} or
#'   \code{summary.stanmvreg} objects is called for its side effect and just returns
#'   its input. The \code{as.data.frame} method for \code{summary.stanmstte}
#'   objects converts the matrix to a data.frame, preserving row and column
#'   names but dropping the \code{print}-related attributes.
#'
#' @seealso \code{\link{prior_summary}} to extract or print a summary of the
#'   priors used for a particular model.
#' @importMethodsFrom rstan summary
summary.stanmstte <- function(object, pars = NULL, regex_pars = NULL,
                            probs = NULL, ..., digits = 1) {

  mstte <- is.stanmstte(object)
  surv <- is.stansurv(object)
  pars <- collect_pars(object, pars, regex_pars)

  if (!used.optimizing(object)) {
    args <- list(object = object$stanfit)
    if (!is.null(probs)){
      args$probs <- probs
    }
    out <- object$stan_summary


    if (is.null(pars) && used.variational(object))
      out <- out[!rownames(out) %in% "log-posterior", , drop = FALSE]
    if (!is.null(pars)) {
      pars <- allow_special_parnames(object, pars)
      out <- out[rownames(out) %in% pars, , drop = FALSE]
    }

    out <- out[!grepl(":_NEW_", rownames(out), fixed = TRUE), , drop = FALSE]
    stats <- colnames(out)
    if ("n_eff" %in% stats)
      out[, "n_eff"] <- round(out[, "n_eff"])
    if ("se_mean" %in% stats) # So people don't confuse se_mean and sd
      colnames(out)[stats %in% "se_mean"] <- "mcse"

  } else { # used optimization
    if (!is.null(probs))
      warning("'probs' ignored if for models fit using optimization.",
              call. = FALSE)
    if (is.null(pars)) {
      famname <- family(object)$family
      mark <- names(object$coefficients)
      if (is.gaussian(famname))
        mark <- c(mark, "sigma")
      if (is.nb(famname))
        mark <- c(mark, "reciprocal_dispersion")
    } else {
      mark <- NA
      if ("alpha" %in% pars)
        mark <- c(mark, "(Intercept)")
      if ("beta" %in% pars)
        mark <- c(mark, setdiff(names(object$coefficients), "(Intercept)"))
      mark <- c(mark, setdiff(pars, c("alpha", "beta")))
      mark <- mark[!is.na(mark)]
    }
    out <- object$stan_summary[mark, , drop=FALSE]
  }

  structure(
    out,
    call          = object$call,
    algorithm     = object$algorithm,
    stan_function = object$stan_function,
    family        = family_plus_link(object),
    formula       = formula(object),
    basehaz       = if (surv){
      basehaz_string(get_basehaz(object))
    } else if (mstte){
      lapply(seq_len(object$n_trans), function(i)
        basehaz_string(get_basehaz(object, i)) )
    }
    else NULL,
    posterior_sample_size = posterior_sample_size(object),
    nobs          = if(mstte) {object$nobs} else {nobs(object)},
    # npreds        = if (is_glm) length(coef(object)) else NULL,
    # ngrps         = if (mer)  ngrps(object)   else NULL,
    nevents       = if (surv){
      object$nevents
    } else if (mstte){
      lapply(seq_len(object$n_trans), function(i)
        object$nevents[[i]])
    } else {
      NULL
    },
    nlcens        = if (surv){
      object$nlcens
    } else if (mstte){
      lapply(seq_len(object$n_trans), function(i)
        object$nlcens[[i]])
    } else {
      NULL
    },
    nrcens        = if (surv){
      object$nrcens
    } else if (mstte){
      lapply(seq_len(object$n_trans), function(i)
        object$nrcens[[i]])
    } else {
      NULL
    },
    nicens        = if (surv){
      object$nicens
    } else if (mstte){
      lapply(seq_len(object$n_trans), function(i)
        object$nicens[[i]])
    } else {
      NULL
    },
    ndelayed      = if (surv){
      object$ndelayed
    } else if (mstte){
      lapply(seq_len(object$n_trans), function(i)
        object$ndelayed[[i]])
    } else {
      NULL
    },
    transition_labels = if(mstte) object$transition_labels else NULL,
    n_trans = if(mstte) object$n_trans else NULL,
    print.digits  = digits,
    priors        = object$prior.info,
    class         = c("summary.stanmstte")
  )
}

#' Print method for summary stanidm
#'
#' @rdname summary.stanmstte
#' @export
#' @method print summary.stanmstte
#'
#' @param x An object of class \code{"summary.stanidm"}.
print.summary.stanmstte <- function(x, digits = max(2, attr(x, "print.digits")),
                                  ...) {
  atts <- attributes(x)
  cat("\nModel Info:\n")
  cat("\n function:       ", atts$stan_function)
  cat("\n algorithm:      ", atts$algorithm)
  cat("\n priors:         ", "see help('prior_summary')")
  cat("\n sample:         ", atts$posterior_sample_size, "(posterior sample size)")

  title <- append_title(atts)

  for(i in seq_len(atts$n_trans)){
    cat("\n ---------------- \n")
    cat(title[[i]])
    cat("\n baseline hazard:", atts$basehaz[[i]] )
    cat("\n formula:        ", formula_string(formula(atts)[[i]]$formula) )
    cat("\n observations:   ", atts$nobs[[i]])
    cat("\n events:         ", atts$nevents[[i]], percent_string(atts$nevents[[i]], atts$nobs[[i]]))
    if (atts$nlcens[[i]] > 0)
      cat("\n left censored:  ", atts$nlcens[[i]], percent_string(atts$nlcens[[i]], atts$nobs[[i]]))
    if (atts$nrcens[[i]] > 0)
      cat("\n right censored: ", atts$nrcens[[i]], percent_string(atts$nrcens[[i]], atts$nobs[[i]]))
    if (atts$nicens[[i]] > 0)
      cat("\n interval cens.: ", atts$nicens[[i]], percent_string(atts$nicens[[i]], atts$nobs[[i]]))
    cat("\n delayed entry:  ", yes_no_string(atts$ndelayed[[i]]) )
  }

  cat("\n\nEstimates:\n")
  sel <- which(colnames(x) %in% c("mcse", "n_eff", "Rhat"))

  if (!length(sel)) {
    .printfr(x, digits)
  } else {
    xtemp <- x[, -sel, drop = FALSE]
    colnames(xtemp) <- paste(" ", colnames(xtemp))
    .printfr(xtemp, digits)
    cat("\nDiagnostics:\n")
    mcse_rhat <- format(round(x[, c("mcse", "Rhat"), drop = FALSE], digits),
                        nsmall = digits)
    n_eff <- format(x[, "n_eff", drop = FALSE], drop0trailing = TRUE)
    print(cbind(mcse_rhat, n_eff), quote = FALSE)
    cat("\nNote: Covariates are arranged in order of transition number, 1->3.\n")
    cat("\nFor each parameter, mcse is the Monte Carlo standard error, ",
        "n_eff is a crude measure of effective sample size, ",
        "and Rhat is the potential scale reduction factor on split chains",
        " (at convergence Rhat=1).\n", sep = '')
  }

  invisible(x)
}



#' @rdname summary.stanmstte
#' @method as.data.frame summary.stanmstte
#' @export
as.data.frame.summary.stanmstte <- function(x, ...) {
  as.data.frame(unclass(x), ...)
}


# --- internal --------------------------------------------
# @param basehaz A list with info about the baseline hazard
.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}

basehaz_string <- function(basehaz, break_and_indent = TRUE) {
  nm <- get_basehaz_name(basehaz)
  switch(nm,
         exp      = "exponential",
         weibull  = "weibull",
         gompertz = "gompertz",
         ms       = "M-splines on hazard scale",
         bs       = "B-splines on log hazard scale",
         piecewise= "piecewise constant on log hazard scale",
         NULL)
}

# @param formula formula object
formula_string <- function(formula, break_and_indent = TRUE) {
  coll <- if (break_and_indent) "--MARK--" else " "
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = coll))
  if (!break_and_indent)
    return(char)
  gsub("--MARK--", "\n\t  ", char, fixed = TRUE)
}

# @param numer,denom The numerator and denominator with which to evaluate a %.
percent_string <- function(numer, denom) {
  val <- round(100 * numer / denom, 1)
  paste0("(", val, "%)")
}

# @param x A logical (or a scalar to be evaluated as a logical).
yes_no_string <- function(x) {
  if (x) "yes" else "no"
}

# get name of aux parameter based on family
.aux_name <- function(object) {
  aux <- character()
  return(aux)
}

.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}


# Allow "alpha", "beta", "varying" as shortcuts
#
# @param object stanreg object
# @param pars result of calling collect_pars(object, pars, regex_pars)
allow_special_parnames <- function(object, pars) {
  pars[pars == "varying"] <- "b"
  pars2 <- NA
  if ("alpha" %in% pars)
    pars2 <- c(pars2, "(Intercept)")
  if ("beta" %in% pars) {
    beta_nms <- names(object$coefficients)
    pars2 <- c(pars2, setdiff(beta_nms, "(Intercept)"))
  }
  if ("b" %in% pars) {
    if (is.mer(object)) {
      pars2 <- c(pars2, b_names(rownames(object$stan_summary), value = TRUE))
      pars[pars == "b"] <- NA
    } else {
      warning("No group-specific parameters. 'varying' ignored.",
              call. = FALSE)
    }
  }
  pars2 <- c(pars2, setdiff(pars, c("alpha", "beta", "varying")))
  pars2[!is.na(pars2)]
}

# Family name with link in parenthesis
# @param x stanreg object
# @param ... Optionally include m to specify which submodel for stanmvreg models
family_plus_link <- function(x, ...) {
  if (is.stansurv(x) | is.stanmstte(x)) {
    return(NULL)
  }
  fam <- family(x, ...)
  if (is.character(fam)) {
    stopifnot(identical(fam, x$method))
    fam <- paste0("ordered [", fam, "]")
  } else if (inherits(x, "betareg")) {
    fam <- paste0("beta [",
                  x$family$link,
                  ", link.phi=",
                  x$family_phi$link,
                  "]")
  } else {
    fam <- paste0(fam$family, " [", fam$link, "]")
  }
  return(fam)
}
