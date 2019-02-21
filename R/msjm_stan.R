#' Bayesian joint longitudinal and multi-state models via Stan
#'
#' \if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
#' Fits a shared parameter joint model for longitudinal and multi-state
#' (e.g. competing risks) data under a Bayesian framework using Stan.
#'
#' @export
#' @param formulaLong A two-sided linear formula object describing both the
#'   fixed-effects and random-effects parts of the longitudinal submodel,
#'   similar in vein to formula specification in the \strong{lme4} package
#'   (see \code{\link[lme4]{glmer}} or the \strong{lme4} vignette for details).
#'   Note however that the double bar (\code{||}) notation is not allowed
#'   when specifying the random-effects parts of the formula, and neither
#'   are nested grouping factors (e.g. \code{(1 | g1/g2))} or
#'   \code{(1 | g1:g2)}, where \code{g1}, \code{g2} are grouping factors.
#'   For a multivariate joint model (i.e. more than one longitudinal marker)
#'   this should be a list of such formula objects, with each element
#'   of the list providing the formula for one of the longitudinal submodels.
#' @param dataLong A data frame containing the variables specified in
#'   \code{formulaLong}. If fitting a multivariate joint model, then this can
#'   be either a single data frame which contains the data for all
#'   longitudinal submodels, or it can be a list of data frames where each
#'   element of the list provides the data for one of the longitudinal
#'   submodels.
#' @param formulaMs A list of two-sided formula object describing the event
#'   submodel(s). The left hand side of the formula should be a \code{Surv()}
#'   object. See \code{\link[survival]{Surv}}. If just providing one formula
#'   a simple one outcome model will be fitted. The architecture of the
#'   multi-state model is represented via the \strong{transition} matrix used to
#'    build the dataframe via \code{ms_prepr}.
#' @param dataMs A data frame containing the variables specified in
#'   \code{formulaMs}. The dataframe should contain a variable trans
#'   for transition as obtained via \code{ms_prepr}.
#' @param time_var A character string specifying the name of the variable
#'   in \code{dataLong} which represents time.
#' @param id_var A character string specifying the name of the variable in
#'   \code{dataLong} which distinguishes between individuals. This can be
#'   left unspecified if there is only one grouping factor (which is assumed
#'   to be the individual). If there is more than one grouping factor (i.e.
#'   clustering beyond the level of the individual) then the \code{id_var}
#'   argument must be specified.
#' @param family The family (and possibly also the link function) for the
#'   longitudinal submodel(s). See \code{\link[lme4]{glmer}} for details.
#'   If fitting a multivariate joint model, then this can optionally be a
#'   list of families, in which case each element of the list specifies the
#'   family for one of the longitudinal submodels.
#' @param assoc A character string or character vector specifying the joint
#'   model association structure. Possible association structures that can
#'   be used include: "etavalue" (the default); "etaslope"; "etaauc";
#'   "muvalue"; "muslope"; "muauc"; "shared_b"; "shared_coef"; or "null".
#'   These are described in the \strong{Details} section below. For a multivariate
#'   joint model, different association structures can optionally be used for
#'   each longitudinal submodel by specifying a list of character
#'   vectors, with each element of the list specifying the desired association
#'   structure for one of the longitudinal submodels. Specifying \code{assoc = NULL}
#'   will fit a joint model with no association structure (equivalent
#'   to fitting separate longitudinal and time-to-event models). It is also
#'   possible to include interaction terms between the association term
#'   ("etavalue", "etaslope", "muvalue", "muslope") and observed data/covariates.
#'   It is also possible, when fitting a multivariate joint model, to include
#'   interaction terms between the association terms ("etavalue" or "muvalue")
#'   corresponding to the different longitudinal outcomes. See the
#'   \strong{Details} section as well as the \strong{Examples} below.
#' @param lag_assoc A non-negative scalar specifying the time lag that should be
#'   used for the association structure. That is, the hazard of the event at
#'   time \emph{t} will be assumed to be associated with the value/slope/auc of
#'   the longitudinal marker at time \emph{t-u}, where \emph{u} is the time lag.
#'   If fitting a multivariate joint model, then a different time lag can be used
#'   for each longitudinal marker by providing a numeric vector of lags, otherwise
#'   if a scalar is provided then the specified time lag will be used for all
#'   longitudinal markers. Note however that only one time lag  can be specified
#'   for linking each longitudinal marker to the
#'   event, and that that time lag will be used for all association structure
#'   types (e.g. \code{"etavalue"}, \code{"etaslope"}, \code{"etaauc"},
#'   \code{"muvalue"}, etc) that are specified for that longitudinal marker in
#'   the \code{assoc} argument.
#' @param grp_assoc Character string specifying the method for combining information
#'   across lower level units clustered within an individual when forming the
#'   association structure. This is only relevant when a grouping factor is
#'   specified in \code{formulaLong} that corresponds to clustering within
#'   individuals. This can be specified as either \code{"sum"}, \code{mean},
#'   \code{"min"} or \code{"max"}. For example, specifying \code{grp_assoc = "sum"}
#'   indicates that the association structure should be based on a summation across
#'   the lower level units clustered within an individual, or specifying
#'   \code{grp_assoc = "mean"}  indicates that the association structure
#'   should be based on the mean (i.e. average) taken across the lower level
#'   units clustered within an individual.
#'   So, for example, specifying \code{assoc = "muvalue"}
#'   and \code{grp_assoc = "sum"} would mean that the log hazard at time
#'   \emph{t} for individual \emph{i} would be linearly related to the sum of
#'   the expected values at time \emph{t} for each of the lower level
#'   units (which may be for example tumor lesions) clustered within that
#'   individual.
#' @param basehaz A character string indicating which baseline hazard to use
#'   for the event submodel. Options are a B-splines approximation estimated
#'   for the log baseline hazard (\code{"bs"}, the default), a Weibull
#'   baseline hazard (\code{"weibull"}, the default), or a piecewise
#'   constant baseline hazard (\code{"piecewise"}). (Note however that there
#'   is currently limited post-estimation functionality available for
#'   models estimated using a piecewise constant baseline hazard).
#' @param basehaz_ops A named list specifying options related to the baseline
#'   hazard. Currently this can include: \cr
#'   \describe{
#'     \item{\code{df}}{A positive integer specifying the degrees of freedom
#'     for the B-splines if \code{basehaz = "bs"}, or the number of
#'     intervals used for the piecewise constant baseline hazard if
#'     \code{basehaz = "piecewise"}. The default is 6.}
#'     \item{\code{knots}}{An optional numeric vector specifying the internal knot
#'     locations for the B-splines if \code{basehaz = "bs"}, or the
#'     internal cut-points for defining intervals of the piecewise constant
#'     baseline hazard if \code{basehaz = "piecewise"}. Knots cannot be
#'     specified if \code{df} is specified. If not specified, then the
#'     default is to use \code{df - 4} knots if \code{basehaz = "bs"},
#'     or \code{df - 1} knots if \code{basehaz = "piecewise"}, which are
#'     placed at equally spaced percentiles of the distribution of
#'     observed event times.}
#'   }
#' @param epsilon The half-width of the central difference used to numerically
#'   calculate the derivate when the \code{"etaslope"} association structure
#'   is used.
#' @param qnodes The number of nodes to use for the Gauss-Kronrod quadrature
#'   that is used to evaluate the cumulative hazard in the likelihood function.
#'   Options are 15 (the default), 11 or 7.
#' @param weights Experimental and should be used with caution. The
#'   user can optionally supply a 2-column data frame containing a set of
#'   'prior weights' to be used in the estimation process. The data frame should
#'   contain two columns: the first containing the IDs for each individual, and
#'   the second containing the corresponding weights. The data frame should only
#'   have one row for each individual; that is, weights should be constant
#'   within individuals.
#' @param init The method for generating the initial values for the MCMC.
#'   The default is \code{"prefit"}, which uses those obtained from
#'   fitting separate longitudinal and time-to-event models prior to
#'   fitting the joint model. The separate longitudinal model is a
#'   (possibly multivariate) generalised linear mixed
#'   model estimated using variational bayes. This is achieved via the
#'   \code{\link{stan_mvmer}} function with \code{algorithm = "meanfield"}.
#'   The separate Cox model is estimated using \code{\link[survival]{coxph}}.
#'   This is achieved
#'   using the and time-to-event models prior
#'   to fitting the joint model. The separate models are estimated using the
#'   \code{\link[lme4]{glmer}} and \code{\link[survival]{coxph}} functions.
#'   This should provide reasonable initial values which should aid the
#'   MCMC sampler. Parameters that cannot be obtained from
#'   fitting separate longitudinal and time-to-event models are initialised
#'   using the "random" method for \code{\link[rstan]{stan}}.
#'   However it is recommended that any final analysis should ideally
#'   be performed with several MCMC chains each initiated from a different
#'   set of initial values; this can be obtained by setting
#'   \code{init = "random"}. In addition, other possibilities for specifying
#'   \code{init} are the same as those described for \code{\link[rstan]{stan}}.
#' @param priorLong,priorEvent,priorEvent_assoc The prior distributions for the
#'   regression coefficients in the longitudinal submodel(s), event submodel,
#'   and the association parameter(s). Can be a call to one of the various functions
#'   provided by \pkg{rstanarm} for specifying priors. The subset of these functions
#'   that can be used for the prior on the coefficients can be grouped into several
#'   "families":
#'
#'   \tabular{ll}{
#'     \strong{Family} \tab \strong{Functions} \cr
#'     \emph{Student t family} \tab \code{normal}, \code{student_t}, \code{cauchy} \cr
#'     \emph{Hierarchical shrinkage family} \tab \code{hs}, \code{hs_plus} \cr
#'     \emph{Laplace family} \tab \code{laplace}, \code{lasso} \cr
#'   }
#'
#'   See the \link[=priors]{priors help page} for details on the families and
#'   how to specify the arguments for all of the functions in the table above.
#'   To omit a prior ---i.e., to use a flat (improper) uniform prior---
#'   \code{prior} can be set to \code{NULL}, although this is rarely a good
#'   idea.
#'
#'   \strong{Note:} Unless \code{QR=TRUE}, if \code{prior} is from the Student t
#'   family or Laplace family, and if the \code{autoscale} argument to the
#'   function used to specify the prior (e.g. \code{\link{normal}}) is left at
#'   its default and recommended value of \code{TRUE}, then the default or
#'   user-specified prior scale(s) may be adjusted internally based on the scales
#'   of the predictors. See the \link[=priors]{priors help page} for details on
#'   the rescaling and the \code{\link{prior_summary}} function for a summary of
#'   the priors used for a particular model.
#' @param priorLong_intercept,priorEvent_intercept The prior distributions
#'   for the intercepts in the longitudinal submodel(s) and event submodel.
#'   Can be a call to \code{normal}, \code{student_t} or
#'   \code{cauchy}. See the \link[=priors]{priors help page} for details on
#'   these functions. To omit a prior on the intercept ---i.e., to use a flat
#'   (improper) uniform prior--- \code{prior_intercept} can be set to
#'   \code{NULL}.
#'
#'   \strong{Note:} The prior distribution for the intercept is set so it
#'   applies to the value when all predictors are centered. Moreover,
#'   note that a prior is only placed on the intercept for the event submodel
#'   when a Weibull baseline hazard has been specified. For the B-splines and
#'   piecewise constant baseline hazards there is not intercept parameter that
#'   is given a prior distribution; an intercept parameter will be shown in
#'   the output for the fitted model, but this just corresponds to the
#'   necessary post-estimation adjustment in the linear predictor due to the
#'   centering of the predictiors in the event submodel.
#'
#' @param priorLong_aux The prior distribution for the "auxiliary" parameters
#'   in the longitudinal submodels (if applicable).
#'   The "auxiliary" parameter refers to a different parameter
#'   depending on the \code{family}. For Gaussian models \code{priorLong_aux}
#'   controls \code{"sigma"}, the error
#'   standard deviation. For negative binomial models \code{priorLong_aux} controls
#'   \code{"reciprocal_dispersion"}, which is similar to the
#'   \code{"size"} parameter of \code{\link[stats]{rnbinom}}:
#'   smaller values of \code{"reciprocal_dispersion"} correspond to
#'   greater dispersion. For gamma models \code{priorLong_aux} sets the prior on
#'   to the \code{"shape"} parameter (see e.g.,
#'   \code{\link[stats]{rgamma}}), and for inverse-Gaussian models it is the
#'   so-called \code{"lambda"} parameter (which is essentially the reciprocal of
#'   a scale parameter). Binomial and Poisson models do not have auxiliary
#'   parameters.
#'
#'   \code{priorLong_aux} can be a call to \code{exponential} to
#'   use an exponential distribution, or \code{normal}, \code{student_t} or
#'   \code{cauchy}, which results in a half-normal, half-t, or half-Cauchy
#'   prior. See \code{\link{priors}} for details on these functions. To omit a
#'   prior ---i.e., to use a flat (improper) uniform prior--- set
#'   \code{priorLong_aux} to \code{NULL}.
#'
#'   If fitting a multivariate joint model, you have the option to
#'   specify a list of prior distributions, however the elements of the list
#'   that correspond to any longitudinal submodel which does not have an
#'   auxiliary parameter will be ignored.
#' @param priorEvent_aux The prior distribution for the "auxiliary" parameters
#'   in the event submodel. The "auxiliary" parameters refers to different
#'   parameters depending on the baseline hazard. For \code{basehaz = "weibull"}
#'   the auxiliary parameter is the Weibull shape parameter. For
#'   \code{basehaz = "bs"} the auxiliary parameters are the coefficients for the
#'   B-spline approximation to the log baseline hazard.
#'   For \code{basehaz = "piecewise"} the auxiliary parameters are the piecewise
#'   estimates of the log baseline hazard.
#' @param prior_covariance Cannot be \code{NULL}; see \code{\link{priors}} for
#'   more information about the prior distributions on covariance matrices.
#'   Note however that the default prior for covariance matrices in
#'   \code{stan_jm} is slightly different to that in \code{\link{stan_glmer}}
#'   (the details of which are described on the \code{\link{priors}} page).
#'
#' @details The \code{stan_jm} function can be used to fit a joint model (also
#'   known as a shared parameter model) for longitudinal and time-to-event data
#'   under a Bayesian framework. The underlying
#'   estimation is carried out using the Bayesian C++ package Stan
#'   (\url{http://mc-stan.org/}). \cr
#'   \cr
#'   The joint model may be univariate (with only one longitudinal submodel) or
#'   multivariate (with more than one longitudinal submodel).
#'   For the longitudinal submodel a (possibly multivariate) generalised linear
#'   mixed model is assumed with any of the \code{\link[stats]{family}} choices
#'   allowed by \code{\link[lme4]{glmer}}. If a multivariate joint model is specified
#'   (by providing a list of formulas in the \code{formulaLong} argument), then
#'   the multivariate longitudinal submodel consists of a multivariate generalized
#'   linear model (GLM) with group-specific terms that are assumed to be correlated
#'   across the different GLM submodels. That is, within
#'   a grouping factor (for example, patient ID) the group-specific terms are
#'   assumed to be correlated across the different GLM submodels. It is
#'   possible to specify a different outcome type (for example a different
#'   family and/or link function) for each of the GLM submodels, by providing
#'   a list of \code{\link[stats]{family}} objects in the \code{family}
#'   argument. Multi-level
#'   clustered data are allowed, and that additional clustering can occur at a
#'   level higher than the individual-level (e.g. patients clustered within
#'   clinics), or at a level lower than the individual-level (e.g. tumor lesions
#'   clustered within patients). If the clustering occurs at a level lower than
#'   the individual, then the user needs to indicate how the lower level
#'   clusters should be handled when forming the association structure between
#'   the longitudinal and event submodels (see the \code{grp_assoc} argument
#'   described above). \cr
#'   \cr
#'   For the event submodel a parametric
#'   proportional hazards model is assumed. The baseline hazard can be estimated
#'   using either a cubic B-splines approximation (\code{basehaz = "bs"}, the
#'   default), a Weibull distribution (\code{basehaz = "weibull"}), or a
#'   piecewise constant baseline hazard (\code{basehaz = "piecewise"}).
#'   If the B-spline or piecewise constant baseline hazards are used,
#'   then the degrees of freedom or the internal knot locations can be
#'   (optionally) specified. If
#'   the degrees of freedom are specified (through the \code{df} argument) then
#'   the knot locations are automatically generated based on the
#'   distribution of the observed event times (not including censoring times).
#'   Otherwise internal knot locations can be specified
#'   directly through the \code{knots} argument. If neither \code{df} or
#'   \code{knots} is specified, then the default is to set \code{df} equal to 6.
#'   It is not possible to specify both \code{df} and \code{knots}. \cr
#'   \cr
#'   Time-varying covariates are allowed in both the
#'   longitudinal and event submodels. These should be specified in the data
#'   in the same way as they normally would when fitting a separate
#'   longitudinal model using \code{\link[lme4]{lmer}} or a separate
#'   time-to-event model using \code{\link[survival]{coxph}}. These time-varying
#'   covariates should be exogenous in nature, otherwise they would perhaps
#'   be better specified as an additional outcome (i.e. by including them as an
#'   additional longitudinal outcome in the joint model). \cr
#'   \cr
#'   Bayesian estimation of the joint model is performed via MCMC. The Bayesian
#'   model includes independent priors on the
#'   regression coefficients for both the longitudinal and event submodels,
#'   including the association parameter(s) (in much the same way as the
#'   regression parameters in \code{\link{stan_glm}}) and
#'   priors on the terms of a decomposition of the covariance matrices of the
#'   group-specific parameters.
#'   See \code{\link{priors}} for more information about the priors distributions
#'   that are available. \cr
#'   \cr
#'   Gauss-Kronrod quadrature is used to numerically evaluate the integral
#'   over the cumulative hazard in the likelihood function for the event submodel.
#'   The accuracy of the numerical approximation can be controlled using the
#'   number of quadrature nodes, specified through the \code{qnodes}
#'   argument. Using a higher number of quadrature nodes will result in a more
#'   accurate approximation.
#'
#'   \subsection{Association structures}{
#'   The association structure for the joint model can be based on any of the
#'   following parameterisations:
#'     \itemize{
#'       \item current value of the linear predictor in the
#'         longitudinal submodel (\code{"etavalue"})
#'       \item first derivative (slope) of the linear predictor in the
#'         longitudinal submodel (\code{"etaslope"})
#'       \item the area under the curve of the linear predictor in the
#'         longitudinal submodel (\code{"etaauc"})
#'       \item current expected value of the longitudinal submodel
#'         (\code{"muvalue"})
#'       \item the area under the curve of the expected value from the
#'         longitudinal submodel (\code{"muauc"})
#'       \item shared individual-level random effects (\code{"shared_b"})
#'       \item shared individual-level random effects which also incorporate
#'         the corresponding fixed effect as well as any corresponding
#'         random effects for clustering levels higher than the individual)
#'         (\code{"shared_coef"})
#'       \item interactions between association terms and observed data/covariates
#'         (\code{"etavalue_data"}, \code{"etaslope_data"}, \code{"muvalue_data"},
#'         \code{"muslope_data"}). These are described further below.
#'       \item interactions between association terms corresponding to different
#'         longitudinal outcomes in a multivariate joint model
#'         (\code{"etavalue_etavalue(#)"}, \code{"etavalue_muvalue(#)"},
#'         \code{"muvalue_etavalue(#)"}, \code{"muvalue_muvalue(#)"}). These
#'         are described further below.
#'       \item no association structure (equivalent to fitting separate
#'         longitudinal and event models) (\code{"null"} or \code{NULL})
#'     }
#'   More than one association structure can be specified, however,
#'   not all possible combinations are allowed.
#'   Note that for the lagged association structures baseline values (time = 0)
#'   are used for the instances
#'   where the time lag results in a time prior to baseline. When using the
#'   \code{"etaauc"} or \code{"muauc"} association structures, the area under
#'   the curve is evaluated using Gauss-Kronrod quadrature with 15 quadrature
#'   nodes. By default, \code{"shared_b"} and \code{"shared_coef"} contribute
#'   all random effects to the association structure; however, a subset of the
#'   random effects can be chosen by specifying their indices between parentheses
#'   as a suffix, for example, \code{"shared_b(1)"} or \code{"shared_b(1:3)"} or
#'   \code{"shared_b(1,2,4)"}, and so on. \cr
#'   \cr
#'   In addition, several association terms (\code{"etavalue"}, \code{"etaslope"},
#'   \code{"muvalue"}, \code{"muslope"}) can be interacted with observed
#'   data/covariates. To do this, use the association term's main handle plus a
#'   suffix of \code{"_data"} then followed by the model matrix formula in
#'   parentheses. For example if we had a variable in our dataset for gender
#'   named \code{sex} then we might want to obtain different estimates for the
#'   association between the current slope of the marker and the risk of the
#'   event for each gender. To do this we would specify
#'   \code{assoc = c("etaslope", "etaslope_data(~ sex)")}. \cr
#'   \cr
#'   It is also possible, when fitting  a multivariate joint model, to include
#'   interaction terms between the association terms themselves (this only
#'   applies for interacting \code{"etavalue"} or \code{"muvalue"}). For example,
#'   if we had a joint model with two longitudinal markers, we could specify
#'   \code{assoc = list(c("etavalue", "etavalue_etavalue(2)"), "etavalue")}.
#'   The first element of list says we want to use the value of the linear
#'   predictor for the first marker, as well as it's interaction with the
#'   value of the linear predictor for the second marker. The second element of
#'   the list says we want to also include the expected value of the second marker
#'   (i.e. as a "main effect"). Therefore, the linear predictor for the event
#'   submodel would include the "main effects" for each marker as well as their
#'   interaction. \cr
#'   \cr
#'   There are additional examples in the \strong{Examples} section below.
#'   }
#'
#' @return A \link[=stanreg-objects]{stanjm} object is returned.
#'
#' @seealso \code{\link{stanreg-objects}}, \code{\link{stanmvreg-methods}},
#'   \code{\link{print.stanmvreg}}, \code{\link{summary.stanmvreg}},
#'   \code{\link{posterior_traj}}, \code{\link{posterior_survfit}},
#'   \code{\link{posterior_predict}}, \code{\link{posterior_interval}},
#'   \code{\link{pp_check}}, \code{\link{ps_check}}, \code{\link{stan_mvmer}}.
#'
#' @examples
#' \donttest{
#' #####
#' # Univariate joint model, with association structure based on the
#' # current value of the linear predictor
#' f1 <- stan_jm(formulaLong = logBili ~ year + (1 | id),
#'               dataLong = pbcLong,
#'               formulaMs = Surv(futimeYears, death) ~ sex + trt,
#'               dataMs = pbcSurv,
#'               time_var = "year",
#'               # this next line is only to keep the example small in size!
#'               chains = 1, cores = 1, seed = 12345, iter = 1000)
#' print(f1)
#' summary(f1)
#'
#' #####
#' # Univariate joint model, with association structure based on the
#' # current value and slope of the linear predictor
#' f2 <- stan_jm(formulaLong = logBili ~ year + (year | id),
#'               dataLong = pbcLong,
#'               formulaMs = Surv(futimeYears, death) ~ sex + trt,
#'               dataMs = pbcSurv,
#'               assoc = c("etavalue", "etaslope"),
#'               time_var = "year",
#'               chains = 1, cores = 1, seed = 12345, iter = 1000)
#' print(f2)
#'
#' #####
#' # Univariate joint model, with association structure based on the
#' # lagged value of the linear predictor, where the lag is 2 time
#' # units (i.e. 2 years in this example)
#' f3 <- stan_jm(formulaLong = logBili ~ year + (1 | id),
#'               dataLong = pbcLong,
#'               formulaMs = Surv(futimeYears, death) ~ sex + trt,
#'               dataMs = pbcSurv,
#'               time_var = "year",
#'               assoc = "etavalue", lag_assoc = 2,
#'               chains = 1, cores = 1, seed = 12345, iter = 1000)
#' print(f3)
#'
#' #####
#' # Univariate joint model, where the association structure includes
#' # interactions with observed data. Here we specify that we want to use
#' # an association structure based on the current value of the linear
#' # predictor from the longitudinal submodel (i.e. "etavalue"), but we
#' # also want to interact this with the treatment covariate (trt) from
#' # pbcLong data frame, so that we can estimate a different association
#' # parameter (i.e. estimated effect of log serum bilirubin on the log
#' # hazard of death) for each treatment group
#' f4 <- stan_jm(formulaLong = logBili ~ year + (1 | id),
#'               dataLong = pbcLong,
#'               formulaMs = Surv(futimeYears, death) ~ sex + trt,
#'               dataMs = pbcSurv,
#'               time_var = "year",
#'               assoc = c("etavalue", "etavalue_data(~ trt)"),
#'               chains = 1, cores = 1, seed = 12345, iter = 1000)
#' print(f4)
#'
#' ######
#' # Multivariate joint model, with association structure based
#' # on the current value and slope of the linear predictor in the
#' # first longitudinal submodel and the area under the marker
#' # trajectory for the second longitudinal submodel
#' mv1 <- stan_jm(
#'         formulaLong = list(
#'           logBili ~ year + (1 | id),
#'           albumin ~ sex + year + (year | id)),
#'         dataLong = pbcLong,
#'         formulaMs = Surv(futimeYears, death) ~ sex + trt,
#'         dataMs = pbcSurv,
#'         assoc = list(c("etavalue", "etaslope"), "etaauc"),
#'         time_var = "year",
#'         chains = 1, cores = 1, seed = 12345, iter = 100)
#' print(mv1)
#'
#' #####
#' # Multivariate joint model, where the association structure is formed by
#' # including the expected value of each longitudinal marker (logBili and
#' # albumin) in the linear predictor of the event submodel, as well as their
#' # interaction effect (i.e. the interaction between the two "etavalue" terms).
#' # Note that whether such an association structure based on a marker by
#' # marker interaction term makes sense will depend on the context of your
#' # application -- here we just show it for demostration purposes).
#' mv2 <- stan_jm(
#'         formulaLong = list(
#'           logBili ~ year + (1 | id),
#'           albumin ~ sex + year + (year | id)),
#'         dataLong = pbcLong,
#'         formulaMs = Surv(futimeYears, death) ~ sex + trt,
#'         dataMs = pbcSurv,
#'         assoc = list(c("etavalue", "etavalue_etavalue(2)"), "etavalue"),
#'         time_var = "year",
#'         chains = 1, cores = 1, seed = 12345, iter = 100)
#'
#' #####
#' # Multivariate joint model, with one bernoulli marker and one
#' # Gaussian marker. We will artificially create the bernoulli
#' # marker by dichotomising log serum bilirubin
#' pbcLong$ybern <- as.integer(pbcLong$logBili >= mean(pbcLong$logBili))
#' mv3 <- stan_jm(
#'         formulaLong = list(
#'           ybern ~ year + (1 | id),
#'           albumin ~ sex + year + (year | id)),
#'         dataLong = pbcLong,
#'         formulaMs = Surv(futimeYears, death) ~ sex + trt,
#'         dataMs = pbcSurv,
#'         family = list(binomial, gaussian),
#'         time_var = "year",
#'         chains = 1, cores = 1, seed = 12345, iter = 1000)
#' }
#'
msjm_stan <- function(formulaLong,
                      dataLong,
                      formulaMs,
                      dataMs,
                      time_var,
                      time_start,
                      transition_labels,
                      n_trans,
                      id_var,
                      family = gaussian,
                      assoc = "etavalue",
                      lag_assoc = 0,
                      grp_assoc,
                      epsilon = 1E-5,
                      basehaz = c("bs", "weibull", "piecewise"),
                      basehaz_ops,
                      qnodes = 15,
                      init = "prefit",
                      weights,
                      priorLong = rstanarm::normal(),
                      priorLong_intercept = rstanarm::normal(),
                      priorLong_aux = rstanarm::cauchy(0, 5),
                      priorMs = rstanarm::normal(),
                      priorMs_intercept = rstanarm::normal(),
                      priorMs_aux = rstanarm::cauchy(),
                      priorMs_assoc = rstanarm::normal(),
                      prior_covariance = rstanarm::lkj(),
                      prior_PD = FALSE,
                      algorithm = c("sampling", "meanfield", "fullrank"),
                      adapt_delta = 0.99, max_treedepth = 10L, QR = FALSE,
                      sparse = FALSE, ... ) {
  
  
  #-----------------------------
  # Pre-processing of arguments
  #-----------------------------
  
  # Set seed if specified
  dots <- list(...)
  if ("seed" %in% names(dots))
    set.seed(dots$seed)
  
  if (!requireNamespace("survival"))
    stop("the 'survival' package must be installed to use this function.")
  
  algorithm <- match.arg(algorithm)
  formulaMs <- validate_arg(formulaMs, type = "formula")
  
  if (missing(basehaz_ops)) basehaz_ops <- NULL
  if (missing(weights))     weights     <- NULL
  if (missing(id_var))      id_var      <- NULL
  if (missing(time_var))    time_var    <- NULL
  if (missing(time_start))  time_start  <- NULL
  if (missing(grp_assoc))   grp_assoc   <- NULL
  
  if (!is.null(weights))
    stop("'weights' are not yet implemented.")
  if (QR)
    stop("'QR' decomposition is not yet implemented.")
  if (sparse)
    stop("'sparse' option is not yet implemented.")
  
  if (is.null(time_var))
    stop("'time_var' must be specified.")
  
  # Formula
  formulaLong <- validate_arg(formulaLong, "formula"); M <- length(formulaLong)
  if (M > 3L)
    stop("'msjm_stan' is currently limited to a maximum of 3 longitudinal outcomes.")
  if(!missing(n_trans) & (length(formulaMs) > 1) ){
    if(n_trans != length(formulaMs))
      stop2("Not matching length of formulaMs and n_trans")
  } else if(missing(n_trans)){
    n_trans     <- length(formulaMs)
  }  else {
    formulaMs <-  maybe_broadcast(formulaMs, n_trans)
  }
  basehaz   <- maybe_broadcast(basehaz, n_trans)
  qnodes    <- maybe_broadcast(qnodes, n_trans)
  # Error if args not supplied together
  supplied_together(formulaLong,  dataLong,  error = TRUE)
  supplied_together(formulaMs, dataMs, error = TRUE)
  
  # Determine whether a joint longitudinal-survival model was specified
  is_msjm <- supplied_together(formulaLong, formulaMs) & (n_trans > 1)
  
  if (is_msjm && is.null(time_start))
    stop("'time_start' variable must be specified.")
  
  # Family
  ok_family_classes <- c("function", "family", "character")
  ok_families <- c("binomial", "gaussian", "Gamma",
                   "inverse.gaussian", "poisson", "neg_binomial_2")
  family <- validate_arg(family, ok_family_classes, validate_length = M)
  family <- lapply(family, validate_famlink, ok_families)
  
  # Assoc
  ok_assoc_classes <- c("NULL", "character")
  assoc <- validate_arg(assoc, ok_assoc_classes, validate_length = M)
  
  # Is priorLong* already a list?
  priorLong           <- broadcast_prior(priorLong, M)
  priorLong_intercept <- broadcast_prior(priorLong_intercept, M)
  priorLong_aux       <- broadcast_prior(priorLong_aux, M)
  
  ## Parse formula
  formulaMs_p <- lapply(formulaMs, function(f) parse_formula(formula = f,
                                                            data = dataMs) )
  # Data
  dataLong  <- validate_arg(dataLong, "data.frame", validate_length = M)
  dataMs <- as.data.frame(dataMs)
  validate_n_trans(dataMs, n_trans)
  
  ## Create data for multi-state
  dataMs <-  lapply(seq_len(n_trans), function(i){
    dataMs[dataMs$trans == i,]
  })
  dataMs <- lapply(seq_len(n_trans), function(i) make_model_data(formulaMs_p[[i]]$tf_form, dataMs[[i]] ) )  # row subsetting etc.
  
  pos <- sapply(dataMs, function(d) nrow(d))
  
  # Observation weights
  has_weights <- !is.null(weights)
  
  # Combine meta information
  stub  <- ifelse(is_msjm, "Long", "y")
  stubs <- paste0(stub, seq(M))
  meta  <- nlist(M,
                 is_msjm,
                 id_var,
                 time_var,
                 time_start,
                 basehaz,
                 basehaz_ops,
                 epsilon,
                 qnodes,
                 stub,
                 stubs,
                 auc_qnodes = 15L)
  
  #--------------------------
  # Longitudinal submodel(s)
  #--------------------------
  
  # info for separate longitudinal submodels
  y_mod <- xapply(formulaLong,
                  dataLong,
                  family,
                  stubs,
                  FUN = handle_y_mod)
  
  # construct single cnms list for all longitudinal submodels
  meta$cnms <- cnms <- get_common_cnms(y_mod, stub = stub)
  
  # construct single list with unique levels for each grouping factor
  meta$flevels <- flevels <- get_common_flevels(y_mod)
  
  # ensure id_var is a valid grouping factor in all submodels
  if (is_msjm) {
    id_var  <- check_id_var (y_mod,   id_var)
    id_list <- check_id_list(y_mod,   id_var)
    weights <- check_weights(weights, id_var)
    meta$id_var  <- id_var
    meta$id_list <- id_list
  }
  
  # observation weights
  y_weights <- lapply(y_mod, handle_weights, weights, id_var)
  
  #----------- Prior distributions -----------#
  
  # valid prior distributions
  ok_dists <- nlist("normal",
                    student_t = "t",
                    "cauchy",
                    "hs",
                    "hs_plus",
                    "laplace",
                    "lasso") # disallow product normal
  ok_dists_int <- c(ok_dists[1:3])
  ok_dists_aux <- c(ok_dists[1:3], exponential = "exponential")
  ok_dists_cov <- c("decov", "lkj")
  
  y_vecs <- fetch(y_mod, "y", "y")     # used in autoscaling, response  vector
  x_mats <- fetch(y_mod, "x", "xtemp") # used in autoscaling, predictor matrix
  
  # note: *_user_prior_*_stuff objects are stored unchanged for constructing
  # prior_summary, while *_prior_*_stuff objects are autoscaled
  
  # priors for longitudinal submodels
  y_links <- fetch(y_mod, "family", "link")
  y_nvars <- fetch(y_mod, "x", "K")
  y_user_prior_stuff <- y_prior_stuff <-
    xapply(priorLong,
           nvars = y_nvars,
           link  = y_links,
           FUN   = handle_glm_prior,
           args  = list(default_scale = 2.5,
                        ok_dists = ok_dists))
  
  y_user_prior_intercept_stuff <- y_prior_intercept_stuff <-
    xapply(priorLong_intercept,
           link  = y_links,
           FUN   = handle_glm_prior,
           args  = list(nvars = 1, default_scale = 10,
                        ok_dists = ok_dists_int))
  
  y_user_prior_aux_stuff <- y_prior_aux_stuff <-
    xapply(priorLong_aux,
           FUN   = handle_glm_prior,
           args  = list(nvars = 1, default_scale = 5, link  = NULL,
                        ok_dists = ok_dists_aux))
  
  b_user_prior_stuff <- b_prior_stuff <-
    handle_cov_prior(prior_covariance, cnms = cnms, ok_dists = ok_dists_cov)
  
  # autoscaling of priors
  y_prior_stuff <-
    xapply(y_prior_stuff,
           response   = y_vecs,
           predictors = x_mats,
           family     = family,
           FUN        = autoscale_prior)
  y_prior_intercept_stuff <-
    xapply(y_prior_intercept_stuff,
           response = y_vecs,
           family   = family,
           FUN      = autoscale_prior)
  y_prior_aux_stuff <-
    xapply(y_prior_aux_stuff,
           response = y_vecs,
           family   = family,
           FUN      = autoscale_prior)
  
  # autoscale priors for sd of random effects (for lkj prior only)
  if (b_prior_stuff$prior_dist_name == "lkj") {
    b_prior_stuff <-
      split_cov_prior(b_prior_stuff, cnms = cnms,
                      submodel_cnms = fetch(y_mod, "z", "group_cnms"))
    b_prior_stuff <-
      xapply(names(cnms), FUN = function(nm) {
        z_mats <- fetch(y_mod, "z", "z", nm)
        xapply(b_prior_stuff[[nm]],
               response   = y_vecs,
               predictors = z_mats,
               family     = family,
               FUN        = autoscale_prior)})
  }
  
  #----------- Data for export to Stan -----------#
  
  standata <- list(
    M           = ai(M),
    has_weights = ai(!all(lapply(weights, is.null))),
    family      = fetch_array(y_mod, "family", "mvmer_family"),
    link        = fetch_array(y_mod, "family", "mvmer_link"),
    weights     = aa(numeric(0)), # not yet implemented
    prior_PD    = ai(prior_PD)
  )
  
  # Dimensions
  standata$has_aux <-
    fetch_array(y_mod, "has_aux", pad_length = 3)
  standata$resp_type <-
    fetch_array(y_mod, "y", "resp_type", pad_length = 3)
  standata$intercept_type <-
    fetch_array(y_mod, "intercept_type", "number", pad_length = 3)
  standata$yNobs <-
    fetch_array(y_mod, "x", "N", pad_length = 3)
  standata$yNeta <-
    fetch_array(y_mod, "x", "N", pad_length = 3) # == Nobs for stan_mvmer models
  standata$yK <-
    fetch_array(y_mod, "x", "K", pad_length = 3)
  
  # Response vectors
  Y_integer <- fetch(y_mod, "y", "integer")
  standata$yInt1 <- if (M > 0) Y_integer[[1]] else aa(integer(0))
  standata$yInt2 <- if (M > 1) Y_integer[[2]] else aa(integer(0))
  standata$yInt3 <- if (M > 2) Y_integer[[3]] else aa(integer(0))
  
  Y_real <- fetch(y_mod, "y", "real")
  standata$yReal1 <- if (M > 0) Y_real[[1]] else aa(double(0))
  standata$yReal2 <- if (M > 1) Y_real[[2]] else aa(double(0))
  standata$yReal3 <- if (M > 2) Y_real[[3]] else aa(double(0))
  
  # Population level design matrices
  X <- fetch(y_mod, "x", "xtemp")
  standata$yX1 <- if (M > 0) X[[1]] else matrix(0,0,0)
  standata$yX2 <- if (M > 1) X[[2]] else matrix(0,0,0)
  standata$yX3 <- if (M > 2) X[[3]] else matrix(0,0,0)
  
  X_bar <- fetch(y_mod, "x", "x_bar")
  standata$yXbar1 <- if (M > 0) aa(X_bar[[1]]) else aa(double(0))
  standata$yXbar2 <- if (M > 1) aa(X_bar[[2]]) else aa(double(0))
  standata$yXbar3 <- if (M > 2) aa(X_bar[[3]]) else aa(double(0))
  
  # Data for group specific terms - group factor 1
  b1_varname <- names(cnms)[[1L]] # name of group factor 1
  b1_nvars <- fetch_(y_mod, "z", "nvars", b1_varname,
                     null_to_zero = TRUE, pad_length = 3)
  b1_ngrps <- fetch_(y_mod, "z", "ngrps", b1_varname)
  if (!n_distinct(b1_ngrps) == 1L)
    stop("The number of groups for the grouping factor '",
         b1_varname, "' should be the same in all submodels.")
  
  standata$bN1 <- b1_ngrps[[1L]] + 1L # add padding for _NEW_ group
  standata$bK1 <- sum(b1_nvars)
  standata$bK1_len <- aa(b1_nvars)
  standata$bK1_idx <- get_idx_array(b1_nvars)
  
  Z1 <- fetch(y_mod, "z", "z", b1_varname)
  Z1 <- lapply(Z1, transpose)
  Z1 <- lapply(Z1, convert_null, "matrix")
  standata$y1_Z1 <- if (M > 0) Z1[[1L]] else matrix(0,0,0)
  standata$y2_Z1 <- if (M > 1) Z1[[2L]] else matrix(0,0,0)
  standata$y3_Z1 <- if (M > 2) Z1[[3L]] else matrix(0,0,0)
  
  Z1_id <- fetch(y_mod, "z", "group_list", b1_varname)
  Z1_id <- lapply(Z1_id, groups)
  Z1_id <- lapply(Z1_id, convert_null, "arrayinteger")
  standata$y1_Z1_id <- if (M > 0) Z1_id[[1L]] else aa(integer(0))
  standata$y2_Z1_id <- if (M > 1) Z1_id[[2L]] else aa(integer(0))
  standata$y3_Z1_id <- if (M > 2) Z1_id[[3L]] else aa(integer(0))
  
  # Data for group specific terms - group factor 2
  if (length(cnms) > 1L) {
    # model has a second grouping factor
    b2_varname <- names(cnms)[[2L]] # name of group factor 2
    b2_nvars <- fetch_(y_mod, "z", "nvars", b2_varname,
                       null_to_zero = TRUE, pad_length = 3)
    b2_ngrps <- fetch_(y_mod, "z", "ngrps", b2_varname)
    if (!n_distinct(b2_ngrps) == 1L)
      stop("The number of groups for the grouping factor '",
           b2_varname, "' should be the same in all submodels.")
    standata$bN2 <- b2_ngrps[[1L]] + 1L # add padding for _NEW_ group
    standata$bK2 <- sum(b2_nvars)
    standata$bK2_len <- aa(b2_nvars)
    standata$bK2_idx <- get_idx_array(b2_nvars)
    
    Z2 <- fetch(y_mod, "z", "z", b2_varname)
    Z2 <- lapply(Z2, transpose)
    Z2 <- lapply(Z2, convert_null, "matrix")
    standata$y1_Z2 <- if (M > 0) Z2[[1L]] else matrix(0,0,0)
    standata$y2_Z2 <- if (M > 1) Z2[[2L]] else matrix(0,0,0)
    standata$y3_Z2 <- if (M > 2) Z2[[3L]] else matrix(0,0,0)
    
    Z2_id <- fetch(y_mod, "z", "group_list", b2_varname)
    Z2_id <- lapply(Z2_id, groups)
    Z2_id <- lapply(Z2_id, convert_null, "arrayinteger")
    standata$y1_Z2_id <- if (M > 0) Z2_id[[1L]] else aa(integer(0))
    standata$y2_Z2_id <- if (M > 1) Z2_id[[2L]] else aa(integer(0))
    standata$y3_Z2_id <- if (M > 2) Z2_id[[3L]] else aa(integer(0))
    
  } else {
    # no second grouping factor
    standata$bN2 <- 0L
    standata$bK2 <- 0L
    standata$bK2_len <- aa(rep(0,3L))
    standata$bK2_idx <- get_idx_array(rep(0,3L))
    standata$y1_Z2 <- matrix(0,0,0)
    standata$y2_Z2 <- matrix(0,0,0)
    standata$y3_Z2 <- matrix(0,0,0)
    standata$y1_Z2_id <- aa(integer(0))
    standata$y2_Z2_id <- aa(integer(0))
    standata$y3_Z2_id <- aa(integer(0))
  }
  
  # Priors
  standata$y_prior_dist_for_intercept <-
    fetch_array(y_prior_intercept_stuff, "prior_dist")
  standata$y_prior_mean_for_intercept <-
    fetch_array(y_prior_intercept_stuff, "prior_mean")
  standata$y_prior_scale_for_intercept <-
    fetch_array(y_prior_intercept_stuff, "prior_scale")
  standata$y_prior_df_for_intercept <-
    fetch_array(y_prior_intercept_stuff, "prior_df")
  
  standata$y_prior_dist_for_aux <-
    fetch_array(y_prior_aux_stuff, "prior_dist")
  standata$y_prior_mean_for_aux <-
    fetch_array(y_prior_aux_stuff, "prior_mean")
  standata$y_prior_scale_for_aux <-
    fetch_array(y_prior_aux_stuff, "prior_scale")
  standata$y_prior_df_for_aux <-
    fetch_array(y_prior_aux_stuff, "prior_df")
  
  standata$y_prior_dist <-
    fetch_array(y_prior_stuff, "prior_dist", pad_length = 3)
  
  prior_mean <- fetch(y_prior_stuff, "prior_mean")
  standata$y_prior_mean1 <- if (M > 0) prior_mean[[1]] else aa(double(0))
  standata$y_prior_mean2 <- if (M > 1) prior_mean[[2]] else aa(double(0))
  standata$y_prior_mean3 <- if (M > 2) prior_mean[[3]] else aa(double(0))
  
  prior_scale <- fetch(y_prior_stuff, "prior_scale")
  standata$y_prior_scale1 <- if (M > 0) aa(prior_scale[[1]]) else aa(double(0))
  standata$y_prior_scale2 <- if (M > 1) aa(prior_scale[[2]]) else aa(double(0))
  standata$y_prior_scale3 <- if (M > 2) aa(prior_scale[[3]]) else aa(double(0))
  
  prior_df <- fetch(y_prior_stuff, "prior_df")
  standata$y_prior_df1 <- if (M > 0) prior_df[[1]] else aa(double(0))
  standata$y_prior_df2 <- if (M > 1) prior_df[[2]] else aa(double(0))
  standata$y_prior_df3 <- if (M > 2) prior_df[[3]] else aa(double(0))
  
  # hs priors only
  standata$y_global_prior_scale <- fetch_array(y_prior_stuff, "global_prior_scale")
  standata$y_global_prior_df    <- fetch_array(y_prior_stuff, "global_prior_df")
  standata$y_slab_df            <- fetch_array(y_prior_stuff, "slab_df")
  standata$y_slab_scale         <- fetch_array(y_prior_stuff, "slab_scale")
  
  # Priors for group specific terms
  standata$t <- length(cnms)
  standata$p <- aa(sapply(cnms, length))
  standata$l <- aa(
    sapply(names(cnms), FUN = function(nm) {
      ngrps <- unique(fetch_(y_mod, "z", "ngrps", nm))
      ngrps + 1L # add padding for _NEW_ group
    }))
  standata$q <- sum(standata$p * standata$l)
  
  if (prior_covariance$dist == "decov") {
    
    # data for decov prior
    standata$prior_dist_for_cov     <- b_prior_stuff$prior_dist
    standata$b_prior_shape          <- b_prior_stuff$prior_shape
    standata$b_prior_scale          <- b_prior_stuff$prior_scale
    standata$b_prior_concentration  <- b_prior_stuff$prior_concentration
    standata$b_prior_regularization <- b_prior_stuff$prior_regularization
    standata$len_concentration      <- length(standata$b_prior_concentration)
    standata$len_regularization     <- length(standata$b_prior_regularization)
    standata$len_theta_L            <- sum(choose(standata$p, 2), standata$p)
    
    # pass empty lkj data
    standata$b1_prior_scale          <- aa(rep(0L, standata$bK1))
    standata$b2_prior_scale          <- aa(rep(0L, standata$bK2))
    standata$b1_prior_df             <- aa(rep(0L, standata$bK1))
    standata$b2_prior_df             <- aa(rep(0L, standata$bK2))
    standata$b1_prior_regularization <- 1.0
    standata$b2_prior_regularization <- 1.0
    
  } else if (prior_covariance$dist == "lkj") {
    
    # data for lkj prior
    b1_prior_stuff          <- b_prior_stuff[[b1_varname]]
    b1_prior_dist           <- fetch_     (b1_prior_stuff, "prior_dist")
    b1_prior_scale          <- fetch_array(b1_prior_stuff, "prior_scale")
    b1_prior_df             <- fetch_array(b1_prior_stuff, "prior_df")
    b1_prior_regularization <- fetch_     (b1_prior_stuff, "prior_regularization")
    
    if (n_distinct(b1_prior_dist) > 1L)
      stop2("Bug found: covariance prior should be the same for all submodels.")
    if (n_distinct(b1_prior_regularization) > 1L)
      stop2("Bug found: prior_regularization should be the same for all submodels.")
    
    standata$prior_dist_for_cov <- unique(b1_prior_dist)
    standata$b1_prior_scale     <- b1_prior_scale
    standata$b1_prior_df        <- b1_prior_df
    standata$b1_prior_regularization <- if (length(b1_prior_regularization))
      unique(b1_prior_regularization) else 1.0
    
    if (standata$bK2 > 0) {
      # model has a second grouping factor
      b2_prior_stuff          <- b_prior_stuff[[b2_varname]]
      b2_prior_scale          <- fetch_array(b2_prior_stuff, "prior_scale")
      b2_prior_df             <- fetch_array(b2_prior_stuff, "prior_df")
      b2_prior_regularization <- fetch_     (b2_prior_stuff, "prior_regularization")
      standata$b2_prior_scale <- b2_prior_scale
      standata$b2_prior_df    <- b2_prior_df
      standata$b2_prior_regularization <- unique(b2_prior_regularization)
    } else {
      # model does not have a second grouping factor
      standata$b2_prior_scale          <- aa(double(0))
      standata$b2_prior_df             <- aa(double(0))
      standata$b2_prior_regularization <- 1.0
    }
    
    # pass empty decov data
    standata$len_theta_L            <- 0L
    standata$len_concentration      <- 0L
    standata$len_regularization     <- 0L
    standata$b_prior_shape          <- aa(rep(0L, standata$t))
    standata$b_prior_scale          <- aa(rep(0L, standata$t))
    standata$b_prior_concentration  <- aa(rep(0L, standata$len_concentration))
    standata$b_prior_regularization <- aa(rep(0L, standata$len_regularization))
  }
  # Names for longitudinal submodel parameters
  y_nms_beta <- uapply(y_mod, get_beta_name_ymod)
  y_nms_int  <- uapply(y_mod, get_int_name_ymod)
  y_nms_aux  <- uapply(y_mod, get_aux_name_ymod)
  y_nms_ppd  <- uapply(y_mod, get_ppd_name)
  
  # Names for group specific coefficients ("b pars")
  b_nms <- get_ranef_name(cnms, flevels)
  
  # Names for Sigma matrix
  y_nms_sigma <- get_Sigma_nms(cnms)
  
  #----------------
  # Multistate submodel(s)
  #----------------
  
  if (is_msjm) { # begin msjm block
    
    # fit separate event submodel
    # ms_mod <-  handle_ms_mod(formula = formulaMs,
    #                          data    = dataMs,
    #                          meta    = meta)
    # Handle time start
    t_start <- lapply(dataMs, function(d) handle_timestart(d, time_start))
    
    ms_mod <- mapply(FUN = handle_e_mod2,
                     t_start = t_start,
                     formula = formulaMs_p,
                     data = dataMs,
                     qnodes = meta$qnodes,
                     basehaz = meta$basehaz,
                     MoreArgs = nlist(meta),
                     SIMPLIFY = FALSE)
    
    meta$has_icens <- lapply(ms_mod, function(ms_mod) ms_mod$has_icens)
    
    meta$id_state <- lapply(dataMs, function(x) extract_id(x, id_var))
    
    # observation weights
    #e_weights <- handle_weights(e_mod, weights, id_var)
    
    assoc_obs <- lapply(seq_len(M), function(m)
      mapply(FUN = assoc_time,
             mod = ms_mod,
             id_state = meta$id_state,
             MoreArgs = nlist(x = dataLong[[m]],
                              t_var = time_var)
             ))
    
    for(j in seq_len(n_trans)){
      ms_mod[[j]]$assoc_obs <- assoc_obs[[1]][[j]]
    }
    
    # check longitudinal observation times are not later than the event time
    for(j in seq_len(n_trans)){
      if(is.data.frame(dataLong)){
        validate_observation_times(
          data = dataLong[assoc_obs[[m]][[j]], ],
          exittime = ms_mod[[j]]$t_end + t_start[[j]],
          id_var     = id_var,
          time_var   = time_var
        )
      } else {
        lapply(seq_len(M), function(m){
          validate_observation_times(
            data = dataLong[[m]][assoc_obs[[m]][[j]], ],
            exittime = ms_mod[[j]]$t_end + t_start[[j]],
            id_var     = id_var,
            time_var   = time_var
          )
        })
      }
    }
    
    #----------- Prior distributions -----------#
    
    # valid prior distributions
    ok_dists_e_aux <- ok_dists[1:3]
    
    # note: *_user_prior_*_stuff objects are stored unchanged for constructing
    # prior_summary, while *_prior_*_stuff objects are autoscaled
    
    # priors for event submodel
    ms_user_prior_stuff <- ms_prior_stuff <- lapply(seq_len(n_trans), function(j){
      handle_glm_prior(priorMs[[j]],
                       nvars = ms_mod[[j]]$K,
                       default_scale = 2.5,
                       link = NULL,
                       ok_dists = ok_dists)
    }
    )
    
    ms_user_prior_intercept_stuff <- ms_prior_intercept_stuff <-
      lapply(priorMs_intercept,
             FUN = handle_glm_prior,
             nvars = 1,
             default_scale = 20,
             link = NULL,
             ok_dists = ok_dists_int)
    
    ms_user_prior_aux_stuff <- ms_prior_aux_stuff <-
      lapply(seq_len(n_trans), function(j) handle_glm_prior(
        priorMs_aux[[j]],
        nvars = ms_mod[[j]]$basehaz$nvars,
        default_scale = get_default_aux_scale(basehaz[[j]]),
        link = NULL,
        ok_dists = ok_dists_e_aux)
      )
    
    # stop null priors if prior_PD is TRUE
    if (prior_PD) {
      if (any(is.null(priorMs)))
        stop("'priorEvent' cannot be NULL if 'prior_PD' is TRUE")
      if (any(is.null(priorMs_intercept)) && ms_mod$has_intercept)
        stop("'priorEvent_intercept' cannot be NULL if 'prior_PD' is TRUE")
      if (any(is.null(prior_aux)))
        stop("'priorEvent_aux' cannot be NULL if 'prior_PD' is TRUE")
    }
    
    # autoscaling of priors
    ms_prior_stuff           <- lapply(seq_len(n_trans), function(k)
      autoscale_prior(ms_prior_stuff[[k]], predictors = ms_mod[[k]]$x) )
    
    ms_prior_intercept_stuff <- lapply(seq_len(n_trans), function(k)
      autoscale_prior(ms_prior_intercept_stuff[[k]]))
    
    ms_prior_aux_stuff       <- lapply(seq_len(n_trans), function(k)
      autoscale_prior(ms_prior_aux_stuff[[k]]) )
    
    #----------- Data for export to Stan -----------#
    # Temporary model with only three states, aka idm
    # Transition 01
    # dimensions
    standata$e_K01              <- ai(ms_mod[[1]]$K)
    standata$nevents01          <- ai(ms_mod[[1]]$nevents)
    #standata$Npat             <- ai(ms_mod[[1]]$Npat)
    #standata$Nevents          <- ai(ms_mod[[1]]$Nevents)
    standata$qnodes01           <- ai(ms_mod[[1]]$qnodes)
    standata$len_cpts01         <- ai(ms_mod[[1]]$len_cpts)
    standata$len_epts01         <- ai(ms_mod[[1]]$len_epts)
    standata$len_qpts01         <- ai(ms_mod[[1]]$len_qpts)
    standata$len_ipts01         <- ai(ms_mod[[1]]$len_ipts)
    standata$idx_cpts01         <- am(ms_mod[[1]]$idx_cpts)
    standata$e_has_intercept01  <- ai(ms_mod[[1]]$has_intercept)
    
    # design matrices
    standata$cpts01             <- aa(ms_mod[[1]]$cpts)
    standata$epts01             <- aa(ms_mod[[1]]$epts)
    standata$qpts01             <- aa(ms_mod[[1]]$qpts)
    standata$ipts01             <- aa(ms_mod[[1]]$ipts)
    standata$qwts01             <- aa(ms_mod[[1]]$qwts)
    standata$iwts01             <- aa(ms_mod[[1]]$iwts)
    standata$eids01             <- aa(ms_mod[[1]]$eids)
    standata$qids01             <- aa(ms_mod[[1]]$qids)
    standata$iids01             <- aa(ms_mod[[1]]$iids)
    standata$e_x01              <- ms_mod[[1]]$x_cpts
    standata$e_xbar01           <- aa(ms_mod[[1]]$x_bar)
    
    # baseline hazard
    standata$basehaz_type01     <- ai(ms_mod[[1]]$basehaz$type)
    standata$basehaz_nvars01    <- ai(ms_mod[[1]]$basehaz$nvars)
    standata$basis_epts01       <- ms_mod[[1]]$basis_epts
    standata$basis_qpts01       <- ms_mod[[1]]$basis_qpts
    standata$basis_ipts01       <- ms_mod[[1]]$basis_ipts
    standata$norm_const01       <- ms_mod[[1]]$norm_const
    
    # priors
    standata$e_prior_dist01              <- ms_prior_stuff[[1]]$prior_dist
    standata$e_prior_dist_for_intercept01 <- ms_prior_intercept_stuff[[1]]$prior_dist
    standata$e_prior_dist_for_aux01      <- ms_prior_aux_stuff[[1]]$prior_dist
    
    # hyperparameters for event submodel priors
    standata$e_prior_mean01               <- ms_prior_stuff[[1]]$prior_mean
    standata$e_prior_scale01              <- ms_prior_stuff[[1]]$prior_scale
    standata$e_prior_df01                 <- ms_prior_stuff[[1]]$prior_df
    standata$e_prior_mean_for_intercept01 <- c(ms_prior_intercept_stuff[[1]]$prior_mean)
    standata$e_prior_scale_for_intercept01 <- c(ms_prior_intercept_stuff[[1]]$prior_scale)
    standata$e_prior_df_for_intercept01   <- c(ms_prior_intercept_stuff[[1]]$prior_df)
    standata$e_prior_mean_for_aux01       <- ms_prior_aux_stuff[[1]]$prior_mean
    standata$e_prior_scale_for_aux01      <- ms_prior_aux_stuff[[1]]$prior_scale
    standata$e_prior_df_for_aux01         <- ms_prior_aux_stuff[[1]]$prior_df
    standata$e_global_prior_scale01       <- ms_prior_stuff[[1]]$global_prior_scale
    standata$e_global_prior_df01          <- ms_prior_stuff[[1]]$global_prior_df
    standata$e_slab_df01                  <- ms_prior_stuff[[1]]$slab_df
    standata$e_slab_scale01               <- ms_prior_stuff[[1]]$slab_scale
    # Transition 02
    # dimensions
    standata$e_K02              <- ai(ms_mod[[2]]$K)
    standata$nevents02          <- ai(ms_mod[[2]]$nevents)
    #standata$Npat             <- ai(ms_mod[[2]]$Npat)
    #standata$Nevents          <- ai(ms_mod[[2]]$Nevents)
    standata$qnodes02           <- ai(ms_mod[[2]]$qnodes)
    standata$len_cpts02         <- ai(ms_mod[[2]]$len_cpts)
    standata$len_epts02         <- ai(ms_mod[[2]]$len_epts)
    standata$len_qpts02         <- ai(ms_mod[[2]]$len_qpts)
    standata$len_ipts02         <- ai(ms_mod[[2]]$len_ipts)
    standata$idx_cpts02         <- am(ms_mod[[2]]$idx_cpts)
    standata$e_has_intercept02  <- ai(ms_mod[[2]]$has_intercept)
    
    # design matrices
    standata$cpts02             <- aa(ms_mod[[2]]$cpts)
    standata$epts02             <- aa(ms_mod[[2]]$epts)
    standata$qpts02             <- aa(ms_mod[[2]]$qpts)
    standata$ipts02             <- aa(ms_mod[[2]]$ipts)
    standata$qwts02             <- aa(ms_mod[[2]]$qwts)
    standata$iwts02             <- aa(ms_mod[[2]]$iwts)
    standata$eids02             <- aa(ms_mod[[2]]$eids)
    standata$qids02             <- aa(ms_mod[[2]]$qids)
    standata$iids02             <- aa(ms_mod[[2]]$iids)
    standata$e_x02              <- ms_mod[[2]]$x_cpts
    standata$e_xbar02           <- aa(ms_mod[[2]]$x_bar)
    
    # baseline hazard
    standata$basehaz_type02     <- ai(ms_mod[[2]]$basehaz$type)
    standata$basehaz_nvars02    <- ai(ms_mod[[2]]$basehaz$nvars)
    standata$basis_epts02       <- ms_mod[[2]]$basis_epts
    standata$basis_qpts02       <- ms_mod[[2]]$basis_qpts
    standata$basis_ipts02       <- ms_mod[[2]]$basis_ipts
    standata$norm_const02       <- ms_mod[[2]]$norm_const
    
    # priors
    standata$e_prior_dist02              <- ms_prior_stuff[[2]]$prior_dist
    standata$e_prior_dist_for_intercept02 <- ms_prior_intercept_stuff[[2]]$prior_dist
    standata$e_prior_dist_for_aux02      <- ms_prior_aux_stuff[[2]]$prior_dist
    
    # hyperparameters for event submodel priors
    standata$e_prior_mean02               <- ms_prior_stuff[[2]]$prior_mean
    standata$e_prior_scale02              <- ms_prior_stuff[[2]]$prior_scale
    standata$e_prior_df02                 <- ms_prior_stuff[[2]]$prior_df
    standata$e_prior_mean_for_intercept02 <- c(ms_prior_intercept_stuff[[2]]$prior_mean)
    standata$e_prior_scale_for_intercept02 <- c(ms_prior_intercept_stuff[[2]]$prior_scale)
    standata$e_prior_df_for_intercept02   <- c(ms_prior_intercept_stuff[[2]]$prior_df)
    standata$e_prior_mean_for_aux02       <- ms_prior_aux_stuff[[2]]$prior_mean
    standata$e_prior_scale_for_aux02      <- ms_prior_aux_stuff[[2]]$prior_scale
    standata$e_prior_df_for_aux02         <- ms_prior_aux_stuff[[2]]$prior_df
    standata$e_global_prior_scale02       <- ms_prior_stuff[[2]]$global_prior_scale
    standata$e_global_prior_df02          <- ms_prior_stuff[[2]]$global_prior_df
    standata$e_slab_df02                  <- ms_prior_stuff[[2]]$slab_df
    standata$e_slab_scale02               <- ms_prior_stuff[[2]]$slab_scale
    # Transition 12
    # dimensions
    standata$e_K12              <- ai(ms_mod[[3]]$K)
    standata$nevents12          <- ai(ms_mod[[3]]$nevents)
    #standata$Npat             <- ai(ms_mod[[3]]$Npat)
    #standata$Nevents          <- ai(ms_mod[[3]]$Nevents)
    standata$qnodes12           <- ai(ms_mod[[3]]$qnodes)
    standata$len_cpts12         <- ai(ms_mod[[3]]$len_cpts)
    standata$len_epts12         <- ai(ms_mod[[3]]$len_epts)
    standata$len_qpts12         <- ai(ms_mod[[3]]$len_qpts)
    standata$len_ipts12         <- ai(ms_mod[[3]]$len_ipts)
    standata$idx_cpts12         <- am(ms_mod[[3]]$idx_cpts)
    standata$e_has_intercept12  <- ai(ms_mod[[3]]$has_intercept)
    
    # design matrices
    standata$cpts12             <- aa(ms_mod[[3]]$cpts)
    standata$epts12             <- aa(ms_mod[[3]]$epts)
    standata$qpts12             <- aa(ms_mod[[3]]$qpts)
    standata$ipts12             <- aa(ms_mod[[3]]$ipts)
    standata$qwts12             <- aa(ms_mod[[3]]$qwts)
    standata$iwts12             <- aa(ms_mod[[3]]$iwts)
    standata$eids12             <- aa(ms_mod[[3]]$eids)
    standata$qids12             <- aa(ms_mod[[3]]$qids)
    standata$iids12             <- aa(ms_mod[[3]]$iids)
    standata$e_x12              <- ms_mod[[3]]$x_cpts
    standata$e_xbar12           <- aa(ms_mod[[3]]$x_bar)
    
    # baseline hazard
    standata$basehaz_type12     <- ai(ms_mod[[3]]$basehaz$type)
    standata$basehaz_nvars12    <- ai(ms_mod[[3]]$basehaz$nvars)
    standata$basis_epts12       <- ms_mod[[3]]$basis_epts
    standata$basis_qpts12       <- ms_mod[[3]]$basis_qpts
    standata$basis_ipts12       <- ms_mod[[3]]$basis_ipts
    standata$norm_const12       <- ms_mod[[3]]$norm_const
    
    # priors
    standata$e_prior_dist12              <- ms_prior_stuff[[3]]$prior_dist
    standata$e_prior_dist_for_intercept12 <- ms_prior_intercept_stuff[[3]]$prior_dist
    standata$e_prior_dist_for_aux12      <- ms_prior_aux_stuff[[3]]$prior_dist
    
    # hyperparameters for event submodel priors
    standata$e_prior_mean12               <- ms_prior_stuff[[3]]$prior_mean
    standata$e_prior_scale12              <- ms_prior_stuff[[3]]$prior_scale
    standata$e_prior_df12                 <- ms_prior_stuff[[3]]$prior_df
    standata$e_prior_mean_for_intercept12 <- c(ms_prior_intercept_stuff[[3]]$prior_mean)
    standata$e_prior_scale_for_intercept12 <- c(ms_prior_intercept_stuff[[3]]$prior_scale)
    standata$e_prior_df_for_intercept12   <- c(ms_prior_intercept_stuff[[3]]$prior_df)
    standata$e_prior_mean_for_aux12       <- ms_prior_aux_stuff[[3]]$prior_mean
    standata$e_prior_scale_for_aux12      <- ms_prior_aux_stuff[[3]]$prior_scale
    standata$e_prior_df_for_aux12         <- ms_prior_aux_stuff[[3]]$prior_df
    standata$e_global_prior_scale12       <- ms_prior_stuff[[3]]$global_prior_scale
    standata$e_global_prior_df12          <- ms_prior_stuff[[3]]$global_prior_df
    standata$e_slab_df12                  <- ms_prior_stuff[[3]]$slab_df
    standata$e_slab_scale12               <- ms_prior_stuff[[3]]$slab_scale
    
    #-----------------------
    # Association structure
    #-----------------------
    
    # define the valid types of association structure
    # !! if order is changed here, then must also change standata$has_assoc !!
    ok_assoc <- c("null",
                  "etavalue",
                  "etaslope",
                  "etaauc",
                  "muvalue",
                  "muslope",
                  "muauc",
                  "shared_b",
                  "shared_coef")
    
    meta$ok_assoc       <- ok_assoc
    meta$ok_assoc_data  <- ok_assoc[c(2,3,5,6)] # ok to interact with covariates
    meta$ok_assoc_int   <- ok_assoc[c(2,5)]     # ok to interact across biomarkers
    meta$ok_assoc_icens <- ok_assoc[c(1:3,5,6)] # ok to use with interval censoring
    
    # check lag time is valid
    lag_assoc <- validate_lag_assoc(lag_assoc, M)
    
    # return an array summarising the association structure
    assoc <- mapply(handle_assoc,
                    user_assoc = assoc,
                    user_lag   = lag_assoc,
                    y_mod      = y_mod,
                    MoreArgs   = nlist(meta))
    assoc <- check_order_of_assoc_interactions(assoc, meta$ok_assoc_int)
    assoc <- set_colnames(assoc, stubs)
    
    # for each submodel, identify grouping factors clustered within 'id_var'
    
    # (i.e. lower level clustering)
    ok_assoc_grp <- c("sum",                      # valid inputs to 'grp_assoc' arg
                      "mean",
                      "min",
                      "max")
    
    ok_assoc_with_grp <- c("etavalue",            # ok to use with non-NULL grp_assoc
                           "etavalue_data",
                           "etaslope",
                           "etaslope_data",
                           "muvalue",
                           "muvalue_data")
    
    grp_basic <- xapply(FUN    = get_basic_grp_info,
                        y_mod  = y_mod,
                        id_var = id_var)
    grp_stuff <- xapply(FUN    = get_extra_grp_info,
                        basic_info = grp_basic,
                        flist  = fetch(y_mod, "z", "group_list"),
                        args   = nlist(id_var, grp_assoc, ok_assoc_grp))
    has_grp <- fetch_(grp_stuff, "has_grp")
    if (not.null(grp_assoc) && !any(has_grp))
      stop2("'grp_assoc' can only be specified when there is a grouping ",
            "factor clustered within patients.")
    if (any(has_grp))
      validate_assoc_with_grp(grp_stuff, assoc, ok_assoc_with_grp)
    
    # design matrices for longitudinal submodel at the quadrature points
    auc_qnodes <- meta$auc_qnodes <- maybe_broadcast(15L, n_trans)
    a_mod <- list()
    
    for(j in (seq_len(n_trans))){
      if(is.data.frame(dataLong)){
        a_mod[[j]] <- lapply(seq_len(M), function(m){
          handle_assocmod2(
            data = dataLong[assoc_obs[[m]][[j]], ],
            assoc     = apply(assoc, 2L, c)[[m]], # converts array to list
            y_mod     = y_mod[[m]],
            grp_stuff = grp_stuff[[m]],
            e_mod = ms_mod[[j]],
            meta = meta,
            j = j)
        })
      } else {
        a_mod[[j]] <- lapply(seq_len(M), function(m){
          handle_assocmod2(
            data = dataLong[[m]][assoc_obs[[m]][[j]], ],
            assoc     = apply(assoc, 2L, c)[[m]], # converts array to list
            y_mod     = y_mod[[m]],
            grp_stuff = grp_stuff[[m]],
            e_mod = ms_mod[[j]],
            meta = meta,
            j = j
          )
        })
      }
    }
    
    # number of association parameters
    a_K <- lapply(a_mod, function(a_mod) get_num_assoc_pars(assoc, a_mod))
    
    # use a stan_mvmer variational bayes model fit for:
    # - obtaining initial values for joint model parameters
    # - obtaining appropriate scaling for priors on association parameters
    dropargs  <- c("chains", "cores", "iter", "refresh", "test_grad", "control")
    init_dots <- list(...); for (i in dropargs) init_dots[[i]] <- NULL
    init_mod  <- stanmodels$mvmer2
    init_data <- standata
    init_pars <- pars_to_monitor2(standata, is_msjm = FALSE)
    init_args <- nlist(object = init_mod,
                       data   = init_data,
                       pars   = init_pars,
                       algorithm = "meanfield")
    init_args[names(init_dots)] <- init_dots
    utils::capture.output(init_fit <- do.call(rstan::vb, init_args))
    init_nms_all <- c(y_nms_int,
                      y_nms_beta,
                      b_nms,
                      y_nms_aux,
                      y_nms_sigma,
                      y_nms_ppd,
                      "log-posterior")
    
    
    
    assoc_ids <- lapply(seq_len(n_trans), function(j)
      mapply( FUN = assoc_id,
              longdata = dataLong,
              MoreArgs = nlist(
                msdata = dataMs[[j]],
                id_var)) )
    
    init_fit  <- replace_stanfit_nms(init_fit, init_nms_all)
    init_mat  <- t(colMeans(am(init_fit))) # posterior means
    init_nms  <- collect_nms(colnames(init_mat), M, stub = "Long")
    init_beta <- lapply(1:M, function(m) init_mat[, init_nms$y[[m]]])
    init_b    <- lapply(seq_len(n_trans), function(j){
      lapply(1:M, function(m) {
        # can drop _NEW_ groups since they are not required for generating
        # the assoc_terms that are used in scaling the priors for
        # the association parameters (ie. the Zt matrix returned by the
        # function 'make_assoc_parts_for_stan' will not be padded).
        b <- init_mat[, init_nms$y_b[[m]]]
        b <- b[!grepl("_NEW_", names(b), fixed = TRUE)]
        b[assoc_ids[[j]]]
      })
    })
    
    if (is.character(init) && (init =="prefit"))
      init <- get_prefit_inits_idm(
        ms_mod = ms_mod,
        assoc_ids_l = assoc_ids,
        transitions = transition_labels,
        prior_stuff_l = ms_prior_stuff,
        prior_aux_stuff_l = ms_prior_aux_stuff,
        standata = standata,
        init_fit = init_fit)
    
    
    
    #----------- Prior distributions -----------#
    
    # Priors for association parameters
    e_user_prior_assoc_stuff <- e_prior_assoc_stuff <-
      mapply( FUN = handle_glm_prior,
              prior = priorEvent_assoc,
              nvars = a_K,
              MoreArgs = nlist(default_scale = 2.5,
                               link          = NULL,
                               ok_dists      = ok_dists),
              SIMPLIFY = FALSE
      )
    
    # Autoscaling of priors
    e_prior_assoc_stuff <- mapply(
      FUN = autoscale_prior2,
      a_k = a_K,
      prior_stuff = e_prior_assoc_stuff,
      parts  = a_mod,
      b      = init_b,
      MoreArgs = nlist(
        family,
        assoc,
        beta   = init_beta
      )
    )
    
    
    #----------- Data for export to Stan -----------#
    
    # dimensions
    standata$assoc01 <- ai(a_K[[1]] > 0L) # any association structure
    standata$a_K01   <- ai(a_K[[1]])      # num association parameters
    
    standata$assoc02 <- ai(a_K[[2]] > 0L) # any association structure
    standata$a_K02   <- ai(a_K[[2]])      # num association parameters
    
    standata$assoc12 <- ai(a_K[[3]] > 0L) # any association structure
    standata$a_K12   <- ai(a_K[[3]])      # num association parameters
    
    # indicators for components required to build association terms
    standata$assoc_uses01 <- make_assoc_component_flags(assoc)
    standata$assoc_uses02 <- make_assoc_component_flags(assoc)
    standata$assoc_uses12 <- make_assoc_component_flags(assoc)
    
    # indicators for each possible type of association structure
    # !! must be careful with corresponding use of indexing in stan code !!
    # !! this is determined by the row ordering of the 'assoc' array     !!
    #   1  = ev
    #   2  = es
    #   3  = ea
    #   4  = mv
    #   5  = ms
    #   6  = ma
    #   7  = shared_b
    #   8  = shared_coef
    #   9  = ev_data
    #   10 = es_data
    #   11 = mv_data
    #   12 = ms_data
    #   13 = evev
    #   14 = evmv
    #   15 = mvev
    #   16 = mvmv
    standata$has_assoc01 <- make_assoc_type_flags(assoc)
    standata$has_assoc02 <- make_assoc_type_flags(assoc)
    standata$has_assoc12 <- make_assoc_type_flags(assoc)
    
    # data for calculating eta, slope, auc in GK quadrature
    standata <- standata_add_assoc_grp01   (standata, a_mod = a_mod[[1]], grp_stuff)
    standata <- standata_add_assoc_grp02   (standata, a_mod = a_mod[[2]], grp_stuff)
    standata <- standata_add_assoc_grp12   (standata, a_mod = a_mod[[3]], grp_stuff)
    
    standata <- standata_add_assoc_xz01    (standata, a_mod = a_mod[[1]], meta = meta, assoc = assoc)
    standata <- standata_add_assoc_xz02    (standata, a_mod = a_mod[[2]], meta = meta, assoc = assoc)
    standata <- standata_add_assoc_xz12    (standata, a_mod = a_mod[[3]], meta = meta, assoc = assoc)
    
    standata <- standata_add_assoc_auc01   (standata, a_mod = a_mod[[1]], meta = meta, j = 1)
    standata <- standata_add_assoc_auc02   (standata, a_mod = a_mod[[2]], meta = meta, j = 2)
    standata <- standata_add_assoc_auc12   (standata, a_mod = a_mod[[3]], meta = meta, j = 3)
    
    standata <- standata_add_assoc_extras01(standata, a_mod = a_mod[[1]], assoc = assoc)
    standata <- standata_add_assoc_extras02(standata, a_mod = a_mod[[2]], assoc = assoc)
    standata <- standata_add_assoc_extras12(standata, a_mod = a_mod[[3]], assoc = assoc)
    
    # hyperparameters for assoc parameter priors
    standata$a_prior_dist01         <- e_prior_assoc_stuff[,1]$prior_dist
    standata$a_prior_mean01          <- e_prior_assoc_stuff[,1]$prior_mean
    standata$a_prior_scale01         <- aa(e_prior_assoc_stuff[,1]$prior_scale)
    standata$a_prior_df01            <- e_prior_assoc_stuff[,1]$prior_df
    standata$a_global_prior_scale01  <- e_prior_assoc_stuff[,1]$global_prior_scale
    standata$a_global_prior_df01     <- e_prior_assoc_stuff[,1]$global_prior_df
    standata$a_slab_df01             <- e_prior_assoc_stuff[,1]$slab_df
    standata$a_slab_scale01          <- e_prior_assoc_stuff[,1]$slab_scale
    
    # centering for association terms
    standata$a_xbar01  <- if (a_K[[1]]) e_prior_assoc_stuff[,1]$a_xbar else numeric(0)
    
    # second
    standata$a_prior_dist02         <- e_prior_assoc_stuff[,2]$prior_dist
    standata$a_prior_mean02          <- e_prior_assoc_stuff[,2]$prior_mean
    standata$a_prior_scale02         <- aa(e_prior_assoc_stuff[,2]$prior_scale)
    standata$a_prior_df02            <- e_prior_assoc_stuff[,2]$prior_df
    standata$a_global_prior_scale02  <- e_prior_assoc_stuff[,2]$global_prior_scale
    standata$a_global_prior_df02     <- e_prior_assoc_stuff[,2]$global_prior_df
    standata$a_slab_df02             <- e_prior_assoc_stuff[,2]$slab_df
    standata$a_slab_scale02          <- e_prior_assoc_stuff[,2]$slab_scale
    
    # centering for association terms
    standata$a_xbar02  <- if (a_K[[2]]) e_prior_assoc_stuff[,2]$a_xbar else numeric(0)
    
    # third
    standata$a_prior_dist12         <- e_prior_assoc_stuff[,3]$prior_dist
    standata$a_prior_mean12          <- e_prior_assoc_stuff[,3]$prior_mean
    standata$a_prior_scale12         <- aa(e_prior_assoc_stuff[,3]$prior_scale)
    standata$a_prior_df12            <- e_prior_assoc_stuff[,3]$prior_df
    standata$a_global_prior_scale12  <- e_prior_assoc_stuff[,3]$global_prior_scale
    standata$a_global_prior_df12     <- e_prior_assoc_stuff[,3]$global_prior_df
    standata$a_slab_df12             <- e_prior_assoc_stuff[,3]$slab_df
    standata$a_slab_scale12          <- e_prior_assoc_stuff[,3]$slab_scale
    
    # centering for association terms
    standata$a_xbar12  <- if (a_K[[3]]) e_prior_assoc_stuff[,3]$a_xbar else numeric(0)
    
  } # end msjm block
  
  #---------------
  # Prior summary
  #---------------
  
  # if (is_jm) {
  #   prior_info <- summarize_jm_prior(
  #     user_priorLong                      = y_user_prior_stuff,
  #     user_priorLong_intercept            = y_user_prior_intercept_stuff,
  #     user_priorLong_aux                  = y_user_prior_aux_stuff,
  #     user_prior_covariance               = prior_covariance,
  #     y_has_intercept                     = fetch_(y_mod, "x", "has_intercept"),
  #     y_has_predictors                    = fetch_(y_mod, "x", "K") > 0,
  #     adjusted_priorLong_scale            = fetch(y_prior_stuff, "prior_scale"),
  #     adjusted_priorLong_intercept_scale  = fetch(y_prior_intercept_stuff, "prior_scale"),
  #     adjusted_priorLong_aux_scale        = fetch(y_prior_aux_stuff, "prior_scale"),
  #     family                              = family,
  #     basehaz                             = e_mod$basehaz,
  #     user_priorEvent                     = e_user_prior_stuff,
  #     user_priorEvent_intercept           = e_user_prior_intercept_stuff,
  #     user_priorEvent_aux                 = e_user_prior_aux_stuff,
  #     user_priorEvent_assoc               = e_user_prior_assoc_stuff,
  #     e_has_intercept                     = standata$e_has_intercept,
  #     e_has_predictors                    = standata$e_K > 0,
  #     has_assoc                           = standata$a_K > 0,
  #     adjusted_priorEvent_scale           = e_prior_stuff$prior_scale,
  #     adjusted_priorEvent_intercept_scale = e_prior_intercept_stuff$prior_scale,
  #     adjusted_priorEvent_aux_scale       = e_prior_aux_stuff$prior_scale,
  #     adjusted_priorEvent_assoc_scale     = e_prior_assoc_stuff$prior_scale)
  # } else {
  #   prior_info <- summarize_jm_prior(
  #     user_priorLong                      = y_user_prior_stuff,
  #     user_priorLong_intercept            = y_user_prior_intercept_stuff,
  #     user_priorLong_aux                  = y_user_prior_aux_stuff,
  #     user_prior_covariance               = prior_covariance,
  #     y_has_intercept                     = fetch_(y_mod, "x", "has_intercept"),
  #     y_has_predictors                    = fetch_(y_mod, "x", "K") > 0,
  #     adjusted_priorLong_scale            = fetch(y_prior_stuff, "prior_scale"),
  #     adjusted_priorLong_intercept_scale  = fetch(y_prior_intercept_stuff, "prior_scale"),
  #     adjusted_priorLong_aux_scale        = fetch(y_prior_aux_stuff, "prior_scale"),
  #     family                              = family)
  # }
  
  #-----------
  # Fit model
  #-----------
  
  # obtain stan model code
  stanfit  <- if (is_msjm) stanmodels$msjm else stanmodels$mvmer2
  
  # specify parameters for stan to monitor
  stanpars <- pars_to_monitor2(standata, is_msjm = is_msjm)
  
  # report type of model to user
  txt1 <- if (M == 1) "uni"   else "multi"
  txt2 <- if (is_msjm)  "joint multi-state model with " else "glmer"
  txt3 <- paste0("Fitting a ", txt1, "variate ", txt2, n_trans, " states.\n\n")
  txt4 <- "Please note the warmup may be much slower than later iterations!\n"
  
  # fit model using stan
  cat(txt3)
  
  if (algorithm == "sampling") { # mcmc
    cat(txt4)
    args <- set_jm_sampling_args(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      init   = init,
      cnms   = cnms,
      user_dots = list(...),
      user_adapt_delta   = adapt_delta,
      user_max_treedepth = max_treedepth,
      show_messages = FALSE)
    stanfit <- do.call(rstan::sampling, args)
    return(stanfit)
  } else { # meanfield or fullrank vb
    args <- nlist(
      object = stanfit,
      data   = standata,
      pars   = stanpars,
      algorithm
    )
    args[names(dots)] <- dots
    stanfit <- do.call(rstan::vb, args)
  }

  if (!isTRUE(check_stanfit(stanfit))) 
    return(standata)
  
  # Sigma values in stanmat
  if (prior_covariance$dist == "decov" && standata$len_theta_L)
    stanfit <- evaluate_Sigma(stanfit, cnms)
  
  if (is_msjm) {
    ms_nms_int <- lapply(seq_len(n_trans), function(i) {
      if(ms_mod[[i]]$has_intercept > 0){
        append_trans(get_int_name_basehaz(basehaz[[i]]), i, transition_labels[i])
      } else {
        NULL
      }
    } )
    ms_nms_beta <- lapply(seq_len(n_trans), function(i) {
      if(ms_mod[[i]]$K > 0){
        append_trans(colnames(ms_mod[[i]]$x), i, transition_labels[i])
      } else {
        NULL
      }
    } )
    
    ms_nms_aux <- lapply(seq_len(n_trans), function(i) {
      if(get_basehaz_name(basehaz[[i]]) != "exp"){
        append_trans(get_aux_name_basehaz(basehaz[[i]]), i, transition_labels[i])
      } else {
        NULL
      }
    } )
    
    ms_nms_assoc <- lapply(a_mod, get_assoc_name, assoc)
    ms_nms_assoc <- lapply(seq_len(n_trans), function(i) {
      append_trans(ms_nms_assoc[[i]], i, transition_labels[i])
    })
  } else {
    ms_nms_int  <- NULL
    ms_nms_int   <- NULL
    ms_nms_aux   <- NULL
    ms_nms_assoc <- NULL
  }
  # define new parameter names
  nms_all <- ulist( c( y_nms_int,
                       y_nms_beta,
                       ms_nms_int,
                       ms_nms_beta,
                       ms_nms_assoc,
                       b_nms,
                       y_nms_aux,
                       ms_nms_aux,
                       y_nms_sigma,
                       y_nms_ppd, 
                       "log-posterior"
                      ) )
  
  # substitute new parameter names into 'stanfit' object
  stanfit <- replace_stanfit_nms(stanfit, nms_all)
  
  # combine elements to add to returned structure
  if (!is_msjm) ms_mod <- a_mod <- assoc <- basehaz <- id_var <- grp_stuff <- NULL
  args <- nlist(.Data = stanfit,
                y_mod, 
                ms_mod,
                a_mod,
                assoc,
                basehaz,
                prior_info = NULL, 
                id_var,
                cnms, 
                flevels,
                grp_stuff)
  
  stanfit <- do.call("structure", rm_null(args, recursive = FALSE))
  
  if (algorithm != "optimizing" && !is(stanfit, "stanfit")) return(stanfit)
  y_mod <- attr(stanfit, "y_mod")
  ms_mod <- attr(stanfit, "ms_mod")
  a_mod <- attr(stanfit, "a_mod")
  cnms  <- attr(stanfit, "cnms")
  flevels <- attr(stanfit, "flevels")
  assoc <- attr(stanfit, "assoc")
  id_var <- attr(stanfit, "id_var")
  basehaz    <- attr(stanfit, "basehaz")
  grp_stuff  <- attr(stanfit, "grp_stuff")
  prior_info <- attr(stanfit, "prior_info")
  stanfit <- drop_attributes(stanfit, "y_mod", "ms_mod", "a_mod", "cnms", 
                             "flevels", "assoc", "id_var", "basehaz", 
                             "grp_stuff", "prior_info")
  
  terms <- c(fetch(y_mod, "terms"), lapply(ms_mod, function(m) terms(m$mod)) ) 
  n_yobs <- fetch_(y_mod, "x", "N")
  n_grps <- sapply(flevels, n_distinct)
  n_subjects <- max( uapply(ms_mod, function(m) m$mod$n ) )
  
  fit <- structure(
    nlist(stanfit, 
               formula = nlist(formulaLong, formulaMs),
               family,
               id_var,
               time_var,
               weights, 
               qnodes, 
               basehaz,
               assoc,
               M,
               n_trans,
               transition_labels,
               cnms, 
               flevels,
               n_grps,
               n_subjects, 
               n_yobs, 
               epsilon,
               algorithm, 
               terms,
               glmod = y_mod,
               msmod = ms_mod, 
               assocmod = a_mod, 
               grp_stuff, 
               dataLong, 
               dataMs,
               prior.info = prior_info,
               stan_function = "stanmsjm", 
               call = match.call(expand.dots = TRUE)),
    class = "stanmsjm")
  
  msjmstan(fit)
}


# ------ Internal
# -------------------

handle_timestart <- function(d, t) {
  return(d[[t]])
}
assoc_id <- function(msdata, longdata, id_var){
  unique(longdata[,id_var]) %in% msdata[,id_var]
}
