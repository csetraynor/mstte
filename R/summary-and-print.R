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
