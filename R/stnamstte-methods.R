#' Methods for stanmstte objects
#'
#' The methods documented on this page are actually some of the least important
#' methods defined for \link[=stanreg-objects]{stanreg} objects. The most
#' important methods are documented separately, each with its own page. Links to
#' those pages are provided in the \strong{See Also} section, below.
#'
#' @name stanmstte-methods
#' @aliases VarCorr fixef ranef ngrps sigma
#' @param ... Ignored, except by the \code{update} method. See
#'   \code{\link{update}}.
#'
#' @details The methods documented on this page are similar to the methods
#'   defined for objects of class 'lm', 'glm', 'glmer', etc. However there are a
#'   few key differences:
#'
#' \describe{
#' \item{\code{residuals}}{
#' Residuals are \emph{always} of type \code{"response"} (not \code{"deviance"}
#' residuals or any other type). However, in the case of \code{\link{stan_polr}}
#' with more than two response categories, the residuals are the difference
#' between the latent utility and its linear predictor.
#' }
#' \item{\code{coef}}{
#' Medians are used for point estimates. See the \emph{Point estimates} section
#' in \code{\link{print.stanreg}} for more details.
#' }
#' \item{\code{se}}{
#' The \code{se} function returns standard errors based on
#' \code{\link{mad}}. See the \emph{Uncertainty estimates} section in
#' \code{\link{print.stanreg}} for more details.
#' }
#' \item{\code{confint}}{
#' For models fit using optimization, confidence intervals are returned via a
#' call to \code{\link[stats]{confint.default}}. If \code{algorithm} is
#' \code{"sampling"}, \code{"meanfield"}, or \code{"fullrank"}, the
#' \code{confint} will throw an error because the
#' \code{\link{posterior_interval}} function should be used to compute Bayesian
#' uncertainty intervals.
#' }
#' }
#'
#' @seealso
#' \itemize{
#'  \item The \code{\link[=print.stanreg]{print}},
#'    \code{\link[=summary.stanreg]{summary}}, and \code{\link{prior_summary}}
#'    methods for stanreg objects for information on the fitted model.
#'  \item \code{\link{launch_shinystan}} to use the ShinyStan GUI to explore a
#'    fitted \pkg{rstanarm} model.
#'  \item The \code{\link[=plot.stanreg]{plot}} method to plot estimates and
#'    diagnostics.
#'  \item The \code{\link{pp_check}} method for graphical posterior predictive
#'    checking.
#'  \item The \code{\link{posterior_predict}} and \code{\link{predictive_error}}
#'    methods for predictions and predictive errors.
#'  \item The \code{\link{posterior_interval}} and \code{\link{predictive_interval}}
#'    methods for uncertainty intervals for model parameters and predictions.
#'  \item The \code{\link[=loo.stanreg]{loo}}, \code{\link{kfold}}, and
#'  \code{\link{log_lik}} methods for leave-one-out or K-fold cross-validation,
#'    model comparison, and computing the log-likelihood of (possibly new) data.
#'  \item The \code{\link[=as.matrix.stanreg]{as.matrix}}, \code{as.data.frame},
#'    and \code{as.array} methods to access posterior draws.
#' }
#'
NULL


#' @rdname stanmstte-methods
#' @export
#' @export fixef
#' @importFrom lme4 fixef
#'
fixef.stanmstte <- function(object, ...) {
  coefs <- object$coefficients
  coefs[b_names(names(coefs), invert = TRUE)]
}

#' @rdname stanmstte-methods
#' @export
#' @export fixef
#' @importFrom lme4 fixef
#'
fixef.stanreg <- function(object, ...) {
  coefs <- object$coefficients
  coefs[b_names(names(coefs), invert = TRUE)]
}
