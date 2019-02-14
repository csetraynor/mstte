# Internal model fitting function for models estimated using
# \code{stan_mvmer} or \code{stan_jm}.
#
# See \code{stan_jm} for a description of the arguments to the
# \code{stan_jm.fit} function call.
#
stan_jm.fit <- function(formulaLong          = NULL,
                        dataLong             = NULL,
                        formulaEvent         = NULL,
                        dataEvent            = NULL,
                        time_var,
                        id_var,
                        family               = gaussian,
                        assoc                = "etavalue",
                        lag_assoc            = 0,
                        grp_assoc,
                        epsilon              = 1E-5,
                        basehaz              = c("bs", "weibull", "piecewise"),
                        basehaz_ops,
                        qnodes               = 15,
                        init                 = "prefit",
                        weights,
                        priorLong            = normal(),
                        priorLong_intercept  = normal(),
                        priorLong_aux        = cauchy(0, 5),
                        priorEvent           = normal(),
                        priorEvent_intercept = normal(),
                        priorEvent_aux       = cauchy(),
                        priorEvent_assoc     = normal(),
                        prior_covariance     = lkj(),
                        prior_PD             = FALSE,
                        algorithm            = c("sampling", "meanfield", "fullrank"),
                        adapt_delta          = NULL,
                        max_treedepth        = 10L,
                        QR                   = FALSE,
                        sparse               = FALSE,
                        ...) {
