library(rstanarm)
library(rstan)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)

stanfile <- "scripts/mstte.stan"

stanfit <-
  rstan::stan(stanfile,
              data = standata,
              chains = 4,
              iter = 500,
              pars = stanpars,
              seed = 1328025050
  )


standata <- mstte_stan(formula = formula,
                     data = ebmt.mstate,
                     basehaz = basehaz,
                     prior           = prior,
                     prior_intercept = prior_intercept,
                     prior_aux       = prior_aux,
                     iter = 2000,
                     chains = 4,
                     control = list(adapt_delta = 0.99)
)
