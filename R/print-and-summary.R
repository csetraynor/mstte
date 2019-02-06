#' Print method for stanidm objects
#' @rdname print.stanidm
#' @export
#' @method print stanidm
print.stanidm <- function(x, digits = 3, ...) {
  ##
  cat("\nAnalysis of independent multi-state risks data \n")
  title <- list("\nState transition 0 -> 1 \n",
                "\nState transition 0 -> 2 \n",
                "\nState transition 1 -> 2 \n")
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
    nms <- grep(paste0("_",i), nms, value = TRUE)
    ppd_nms <- grep("^mean_PPD", nms, value = TRUE)
    nms <- setdiff(nms, ppd_nms)
    coef_mat <- mat[, nms, drop = FALSE]
    estimates <- .median_and_madsd(coef_mat)
    if(has_intercept(x$basehaz[[i]])){
      nms_int  <- paste_i( get_int_name_basehaz(get_basehaz(x)[[i]]), i)
    } else {
      nms_int = NULL
    }

    if(get_basehaz_name(x$basehaz[[i]]) != "exp"){
      nms_aux  <- paste_i(get_aux_name_basehaz(get_basehaz(x)[[i]]), i)
    } else {
      nms_aux = NULL
    }
    if(ncol(get_x(x)[[i]]) > 0){
      nms_beta <- paste_i(
        setdiff(gsub(paste0("_.") , "", rownames(estimates)),
                gsub(paste0("_.") , "", c(nms_int, nms_aux)) ), i)
    } else {
      nms_beta <- NULL
    }

    estimates <- suppressWarnings(
      cbind(estimates,
            "exp(Median)" = c(rep(NA, length(nms_int)), exp(estimates[nms_beta, "Median"]), rep(NA, length(nms_aux)))
      )
    )
    .printfr(estimates, digits, ...)

  }

  cat("\n ---------------- \n")
  cat("* For help interpreting the printed output see ?print.stanidm\n")
  cat("* For info on the priors used see ?prior_summary.stanreg\n")
}


#' Summary method for stanidm objects
#'
#' Summaries of parameter estimates and MCMC convergence diagnostics
#' (Monte Carlo error, effective sample size, Rhat).
#' @export
#' @method summary stanidm
summary.stanidm <- function(object, pars = NULL, regex_pars = NULL,
                            probs = NULL, ..., digits = 1) {
  surv <- is.surv(object)
  # mer  <- is.mer(object)
  idm <- is.idm(object)
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

  is_glm <-
    isTRUE(object$stan_function %in% c("stan_glm", "stan_glm.nb", "stan_lm"))

  structure(
    out,
    call          = object$call,
    algorithm     = object$algorithm,
    stan_function = object$stan_function,
    family        = family_plus_link(object),
    formula       = formula(object),
    basehaz       = if (surv){
      basehaz_string(get_basehaz(object))
    } else if (idm){
      lapply(seq_along(object$basehaz), function(i)
        basehaz_string(get_basehaz(object)[[i]]) )
    }
    else NULL,
    posterior_sample_size = posterior_sample_size(object),
    nobs          = if(idm) {object$nobs} else {nobs(object)},
    npreds        = if (is_glm) length(coef(object)) else NULL,
    # ngrps         = if (mer)  ngrps(object)   else NULL,
    nevents       = if (surv){
      object$nevents
    } else if (idm){
      lapply(seq_along(object$basehaz), function(i)
        object$nevents[[i]])
    } else {
      NULL
    },
    nlcens        = if (surv){
      object$nlcens
    } else if (idm){
      lapply(seq_along(object$basehaz), function(i)
        object$nlcens[[i]])
    } else {
      NULL
    },
    nrcens        = if (surv){
      object$nrcens
    } else if (idm){
      lapply(seq_along(object$basehaz), function(i)
        object$nrcens[[i]])
    } else {
      NULL
    },
    nicens        = if (surv){
      object$nicens
    } else if (idm){
      lapply(seq_along(object$basehaz), function(i)
        object$nicens[[i]])
    } else {
      NULL
    },
    ndelayed      = if (surv){
      object$ndelayed
    } else if (idm){
      lapply(seq_along(object$basehaz), function(i)
        object$ndelayed[[i]])
    } else {
      NULL
    },
    print.digits  = digits,
    priors        = object$prior.info,
    class         = c("summary.stanidm")
  )
}

#' Print method for summary stanidm
#'
#' @rdname summary.stanidm
#' @export
#' @method print summary.idm
#'
#' @param x An object of class \code{"summary.stanidm"}.
print.summary.stanidm <- function(x, digits = max(2, attr(x, "print.digits")),
                                  ...) {
  atts <- attributes(x)
  cat("\nModel Info:\n")
  cat("\n function:       ", atts$stan_function)
  cat("\n algorithm:      ", atts$algorithm)
  cat("\n priors:         ", "see help('prior_summary')")
  cat("\n sample:         ", atts$posterior_sample_size, "(posterior sample size)")

  title <- list("\nState transition 1 (0 -> 1) \n",
                "\nState transition 2 (0 -> 2) \n",
                "\nState transition 3 (1 -> 2) \n")
  for(i in seq_along(atts$basehaz)){
    cat("\n ---------------- \n")
    cat(title[[i]])
    cat("\n baseline hazard:", basehaz_string(atts$basehaz[[i]]) )
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



#' @rdname summary.stanidm
#' @method as.data.frame summary.stanidm
#' @export
as.data.frame.summary.stanidm <- function(x, ...) {
  as.data.frame(unclass(x), ...)
}


# internal ----------------------------------------------------------------
.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
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
    beta_nms <- if (is.mer(object))
      names(fixef(object)) else names(object$coefficients)
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
  if (is.stansurv(x) | is.idm(x)) {
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

# @param formula formula object
formula_string <- function(formula, break_and_indent = TRUE) {
  coll <- if (break_and_indent) "--MARK--" else " "
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = coll))
  if (!break_and_indent)
    return(char)
  gsub("--MARK--", "\n\t  ", char, fixed = TRUE)
}

# get name of aux parameter based on family
.aux_name <- function(object) {
  aux <- character()
  if (!is_polr(object) && !is.stansurv(object) && !is.idm(object)) {
    aux <- .rename_aux(family(object))
    if (is.na(aux)) {
      aux <- character()
    }
  }
  return(aux)
}

# print anova table for stan_aov models
# @param x stanreg object created by stan_aov()
print_anova_table <- function(x, digits, ...) {
  labels <- attributes(x$terms)$term.labels
  patterns <- gsub(":", ".*:", labels)
  dnms <- dimnames(extract(x$stanfit, pars = "beta",
                           permuted = FALSE))$parameters
  groups <- sapply(patterns, simplify = FALSE, FUN = grep, x = dnms)
  names(groups) <- gsub(".*", "", names(groups), fixed = TRUE)
  groups <- groups[sapply(groups, length) > 0]
  effects_dim <- dim(x$effects)
  effects <- x$effects^2
  effects <- sapply(groups, FUN = function(i) {
    apply(effects[, , i, drop = FALSE], 1:2, mean)
  })
  dim(effects) <- c(effects_dim[-3], ncol(effects))
  dim(effects) <- c(nrow(effects) * ncol(effects), dim(effects)[3])
  colnames(effects) <- paste("Mean Sq", names(groups))
  anova_table <- .median_and_madsd(effects)
  cat("\nANOVA-like table:\n")
  .printfr(anova_table, digits, ...)
}

# @param basehaz A list with info about the baseline hazard
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

# @param x A logical (or a scalar to be evaluated as a logical).
yes_no_string <- function(x) {
  if (x) "yes" else "no"
}

# @param numer,denom The numerator and denominator with which to evaluate a %.
percent_string <- function(numer, denom) {
  val <- round(100 * numer / denom, 1)
  paste0("(", val, "%)")
}

# equivalent to isFALSE(object$compute_mean_PPD)
no_mean_PPD <- function(object) {
  x <- object$compute_mean_PPD
  is.logical(x) && length(x) == 1L && !is.na(x) && !x
}
