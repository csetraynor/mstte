roxygen2::roxygenise(clean = TRUE)
source("~/aeim/R/simulation.R")

betas01_t = c(trt = -0.2)
betas02_t = c(trt =-0.4)
betas12_t = c(trt =-0.5)
lambdas01_t = 0.325
lambdas02_t = 0.36
lambdas12_t = 0.39
gammas01_t = 1.9
gammas02_t = 2.2
gammas12_t = 2.5
cens = c(4.5, 5.5)

set.seed(9911)
covs <- data.frame(id = 1:100, trt = stats::rbinom(100, 1L, 0.5))

sim_wei <- rsimid(
  dist01 = "weibull",
  dist02 = "weibull",
  dist12 = "weibull",
  betas01 = betas01_t,
  betas02 = betas02_t,
  betas12 = betas12_t,
  lambdas01 = lambdas01_t,
  lambdas02 = lambdas02_t,
  lambdas12 = lambdas12_t,
  gammas01 = gammas01_t,
  gammas02 = gammas02_t,
  gammas12 = gammas12_t,
  x = covs,
  cens = cens
)

#sim_wei$time_diff = sim_wei$os_time - sim_wei$df_time
tmat_mst <- mstate::trans.illdeath(names=c("diagnosis","relapse","death"))
sim_wei_mstate <- mstate::msprep(time=c(NA,"df_time","os_time"),
                         status=c(NA,"df_event","os_event"),
                         data = sim_wei,
                         trans=tmat_mst)
sim_wei_mstate <- dplyr::left_join(sim_wei_mstate,
                                   sim_wei[ , c("id", "trt")])


formula = lapply(1:3, function (x)
  as.formula(Surv(time=time,event=status) ~ trt) )

basehaz = lapply(1:3, function(x)
  "weibull")

prior_intercept = lapply(1:3, function(x)
  rstanarm::normal() )

prior_aux = lapply(1:3, function(x)
  rstanarm::cauchy() )

options(mc.cores = 4)
stanfit <- mstte_stan(formula = formula,
                      data = sim_wei_mstate,
                      transition_labels = c("DP", "DX", "DPDX"),
                      basehaz = basehaz,
                      prior           = lapply(1:3, function(x)
                        rstanarm::normal() ),
                      prior_intercept = prior_intercept,
                      prior_aux       = prior_aux,
                      iter = 1,
                      chains = 1,
                      control = list(adapt_delta = 0.99)
)

saveRDS(nlist(stanfit,sim_wei_mstate ), "~/rfactory/mstte-data/baseline_sim.RDS")
##########################################



library(survival)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
stanfit1 <- idm_stan(formula01 = Surv(time=df_time,event=df_event)~trt,
                     formula02 = Surv(time=os_time,event=os_event)~trt,
                     formula12 = Surv(time=time_diff,event=os_event)~trt,
                     data = sim_wei,
                     basehaz01 = "weibull",
                     basehaz02 = "weibull",
                     basehaz12 = "weibull",
                     prior01           = rstanarm::normal(),
                     prior_intercept01 = rstanarm::normal(),
                     prior_aux01       = rstanarm::normal(),
                     prior02           = rstanarm::normal(),
                     prior_intercept02 = rstanarm::normal(),
                     prior_aux02       = rstanarm::normal(),
                     prior12           = rstanarm::normal(),
                     prior_intercept12 = rstanarm::normal(),
                     prior_aux12       = rstanarm::normal(),
                     iter = 8000,
                     chains = 4,
                     control = list(adapt_delta = 0.99, max_treedepth = 15)
)

exp_out <- as.data.frame(summary(stanfit1))
exp_out <- exp_out[1:6, ]

actual <- c(
  "(Intercept)_1"  = log(lambdas01_t),
  "weibull-shape_1" = gammas01_t,

  "(Intercept)_2" = log(lambdas02_t),
  "weibull-shape_2" = gammas02_t,

  "(Intercept)_3" = log(lambdas12_t),
  "weibull-shape_3" = gammas12_t)

exp_out <- cbind(actual, exp_out)

write.csv(exp_out, "actual_weibull.csv")




#### Bivariate

betas01_t = c(trt = -0.23, age = -0.1)
betas02_t = c(trt =-0.19, age = 0.22)
betas12_t = c(trt =-0.26, age = 0.26)
lambdas01_t = 0.22
lambdas02_t = 0.26
lambdas12_t = 0.29
gammas01_t = 1.9
gammas02_t = 2.2
gammas12_t = 2.5
cens = c(7.5, 9.5)



set.seed(9911)
covs <- data.frame(id = 1:16000,
                   trt = stats::rbino00, 1L, 0.5),
age = rnor400000))

sim_wei <- rsimid(
  dist01 = "weibull",
  dist02 = "weibull",
  dist12 = "weibull",
  betas01 = betas01_t,
  betas02 = betas02_t,
  betas12 = betas12_t,
  lambdas01 = lambdas01_t,
  lambdas02 = lambdas02_t,
  lambdas12 = lambdas12_t,
  gammas01 = gammas01_t,
  gammas02 = gammas02_t,
  gammas12 = gammas12_t,
  x = covs,
  cens = cens
)

sim_wei$time_diff = sim_wei$os_time - sim_wei$df_time

library(survival)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
stanfit2 <- idm_stan(formula01 = Surv(time=df_time,event=df_event)~trt+age,
                     formula02 = Surv(time=os_time,event=os_event)~trt+age,
                     formula12 = Surv(time=time_diff,event=os_event)~trt+age,
                     data = sim_wei,
                     basehaz01 = "weibull",
                     basehaz02 = "weibull",
                     basehaz12 = "weibull",
                     prior01           = rstanarm::normal(),
                     prior_intercept01 = rstanarm::normal(),
                     prior_aux01       = rstanarm::normal(),
                     prior02           = rstanarm::normal(),
                     prior_intercept02 = rstanarm::normal(),
                     prior_aux02       = rstanarm::normal(),
                     prior12           = rstanarm::normal(),
                     prior_intercept12 = rstanarm::normal(),
                     prior_aux12       = rstanarm::normal(),
                     iter = 8000,
                     chains = 4,
                     control = list(adapt_delta = 0.99, max_treedepth = 15)
)





## ----plot-alpha-vs-test-------------------------------------------------- BAYES PLOT -
library(bayesplot)
posterior <- as.array(stanfit)
dim(posterior)
dimnames(posterior)

color_scheme_set("gray")
p1 <- mcmc_hex(posterior, pars = c("(Intercept)_1", "weibull-shape_1")) +
  geom_point(aes(x = log(lambdas01_t), y = gammas01_t ), shape = 25, size = 3, fill = "black") + ggtitle('Transition 0 -> 1')  +
  labs(x = "scale", y = "shape")

p2 <- mcmc_hex(posterior, pars = c("(Intercept)_2", "weibull-shape_2")) +
  geom_point(aes(x = log(lambdas02_t), y = gammas02_t ), shape = 25, size = 3, fill = "black") + ggtitle('Transition 0 -> 2')  +
  labs(x = "scale", y = "shape")

p3 <- mcmc_hex(posterior, pars = c("(Intercept)_3", "weibull-shape_3")) +
  geom_point(aes(x = log(lambdas12_t), y = gammas12_t ), shape = 25, size = 3, fill = "black") + ggtitle('Transition 1 -> 2')  +
  labs(x = "scale", y = "shape")


library(cowplot)
title <- ggdraw() + draw_label("Posterior joint distribution of Weibull-shape and scale showing true parameter values", fontface='bold')

p <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
pp <- cowplot::plot_grid(title, p, ncol = 1, rel_heights=c(0.1, 1))
pp

ggsave(plot = pp, filename = "nipp_wei.png", height = 4, width = 12, units = "in", dpi = 600)


# --------------------------------------------------

exp_out <- as.data.frame(summary(stanfit))
exp_out <- exp_out[1:6, ]

actual <- c(
  "(Intercept)_1"  = log(lambdas01_t),
  "weibull-shape_1" = gammas01_t,

  "(Intercept)_2" = log(lambdas02_t),
  "weibull-shape_2" = gammas02_t,

  "(Intercept)_3" = log(lambdas12_t),
  "weibull-shape_3" = gammas12_t)

exp_out <- cbind(actual, exp_out)

write.csv(exp_out, "actual_weibull.csv")

pp <- pairs(stanfit$stanfit)

ggsave(plot = pp, filename = "pp_wei.png", height = 3.8, width = 7, units = "in", dpi = 600)

lambdas01_t = 0.08
lambdas02_t = 0118
lambdas12_t = 0513
gammas01_t = 0.22
gammas02_t = 0.34
gammas12_t = 0.56
cens = c(4.5, 5.5)

set.seed(9911)
covs <- data.frame(id = 1:20000, trt = stats::rbinom(20000, 1L, 0.5))

sim_gomp <- rsimid(
  dist01 = "gompertz",
  dist02 = "gompertz",
  dist12 = "gompertz",
  lambdas01 = lambdas01_t,
  lambdas02 = lambdas02_t,
  lambdas12 = lambdas12_t,
  gammas01 = gammas01_t,
  gammas02 = gammas02_t,
  gammas12 = gammas12_t,
  x = covs,
  cens = cens
)

sim_gomp$time_diff = sim_gomp$os_time - sim_gomp$df_time

stanfit2 <- idm_stan(formula01 = Surv(time=df_time,event=df_event)~1,
                     formula02 = Surv(time=os_time,event=os_event)~1,
                     formula12 = Surv(time=time_diff,event=os_event)~1,
                     data = sim_gomp,
                     basehaz01 = "gompertz",
                     basehaz02 = "gompertz",
                     basehaz12 = "gompertz",
                     prior01           = rstanarm::normal(),
                     prior_intercept01 = rstanarm::normal(),
                     prior_aux01       = rstanarm::normal(),
                     prior02           = rstanarm::normal(),
                     prior_intercept02 = rstanarm::normal(),
                     prior_aux02       = rstanarm::normal(),
                     prior12           = rstanarm::normal(),
                     prior_intercept12 = rstanarm::normal(),
                     prior_aux12       = rstanarm::normal(),
                     iter = 4000,
                     chains = 4,
                     control = list(adapt_delta = 0.99)
)

gomp_out <- as.data.frame(summary(stanfit2))

gomp_out <- gomp_out[1:6, ]

actual <- c(
  "(Intercept)_1"  = log(lambdas01_t),
  "gompertz-scale_1" = gammas01_t,

  "(Intercept)_2" = log(lambdas02_t),
  "gompertz-scale_2" = gammas02_t,

  "(Intercept)_3" = log(lambdas12_t),
  "gompertz-scale_3" = gammas12_t)

gomp_out <- cbind(actual, gomp_out)

write.csv(gomp_out, "actual_gompertz.csv")

library(cowplot)

## ----sim-extract-alpha---------------------------------------------------
pp_lambda01 <- exp( rstan::extract(stanfit2$stanfit,'alpha01')$alpha01 )
pp_lambda02 <- exp( rstan::extract(stanfit2$stanfit,'alpha02')$alpha02 )
pp_lambda12 <- exp( rstan::extract(stanfit2$stanfit,'alpha12')$alpha12 )


pp_gamma01 <- rstan::extract(stanfit2$stanfit,'aux01')$aux01
pp_gamma02 <- rstan::extract(stanfit2$stanfit,'aux02')$aux02
pp_gamma12 <- rstan::extract(stanfit2$stanfit,'aux12')$aux12

pp_gamma01 <- rstan::extract(fit,'aux01')$aux01
pp_beta01 <- rstan::extract(fit,'beta01')$beta01

## ----plot-alpha-vs-test--------------------------------------------------
library(bayesplot)
posterior <- as.array(stanfit)
dim(posterior)
dimnames(posterior)

color_scheme_set("gray")
mcmc_hex(posterior, pars = c("(Intercept)_1", "weibull-shape_1")) +
  geom_point(aes(x= ))




p1 <- ggplot(data.frame(lambda = pp_lambda01)) +
  geom_density(aes(x = lambda)) +
  geom_vline(aes(xintercept = lambdas01_t), colour = 'red') +
  ggtitle('Posterior distribution of lambda 0 -> 1\nshowing true value')

p2 <- ggplot(data.frame(lambda = pp_lambda02)) +
  geom_density(aes(x = lambda)) +
  geom_vline(aes(xintercept = lambdas02_t), colour = 'red') +
  ggtitle('Posterior distribution of lambda 0 -> 2\nshowing true value')

p3 <- ggplot(data.frame(lambda = pp_lambda12)) +
  geom_density(aes(x = lambda)) +
  geom_vline(aes(xintercept = lambdas12_t), colour = 'red') +
  ggtitle('Posterior distribution of lambda 1 -> 2\nshowing true value')

cowplot::plot_grid(p1, p2, p3)


## ----plot-mu-vs-test-----------------------------------------------------
ggplot(data.frame(lambda = pp_lambda01, gamma = pp_gamma01)) +
  geom_density(aes(x = gamma)) +
  geom_vline(aes(xintercept = gammas01_t), colour = 'red') +
  ggtitle('Posterior distribution of mu\nshowing true value in red')

ggplot(data.frame(lambda = pp_lambda02, gamma = pp_gamma02)) +
  geom_density(aes(x = gamma)) +
  geom_vline(aes(xintercept = gammas02_t), colour = 'red') +
  ggtitle('Posterior distribution of mu\nshowing true value in red')


## ----plot-mu-vs-alpha----------------------------------------------------

pp_joint01 <- pp_lambda01 *  pp_gamma01

library(ggsubplot)
plot_data <- data.frame(lambda = pp_lambda01,
                        gamma = pp_gamma01,
                        joint = pp_joint01)

p1 <- ggplot(plot_data, aes(x=lambda, y=gamma, z = joint)) +
  geom_contour(aes(x = lambda, y = gamma, colour = stat(level)),
               linetype = "dashed") +
  geom_point(aes(x = lambdas01_t, y = gammas01_t),
             colour = 'red',
             size = 2) +
  ggtitle('Transition 0 -> 1')  +
  labs(x = "rate", y = "shape")

#ggtitle('Posterior distributions of rate and shape \nshowing true parameter values')

p2 <- ggplot(data.frame(lambda = pp_lambda02,
                        gamma = pp_gamma02,
                        joint = pp_lambda02 * pp_gamma02),
             aes(x = lambda, y = gamma, z = joint)) +
  geom_contour(aes(x = lambda, y = gamma, colour = stat(level)), linetype = "dashed", colour = "white", n = 100, contour = TRUE) +
  geom_point(aes(x = lambdas02_t, y = gammas02_t), colour = 'red', size = 2) + ggtitle('Transition 0 -> 2')  +
  labs(x = "rate", y = "shape")


# ggtitle('Posterior distributions of lambda and gamma\nshowing true parameter values') +
# labs(x = "rate", y = "sahpe")

p3 <- ggplot(data.frame(lambda = pp_lambda12,
                        gamma = pp_gamma12,
                        joint = pp_lambda12 * pp_gamma12),
             aes(x = lambda, y = gamma, z = joint) ) +
  geom_contour(aes(x = lambda, y = gamma, colour = stat(level)), linetype = "dashed") +
  geom_point(aes(x = lambdas12_t, y = gammas12_t), colour = 'red', size = 2) + ggtitle('Transition 1 -> 2')  +
  labs(x = "rate", y = "shape")
# ggtitle('Posterior distributions of lambda and gamma\nshowing true parameter values')
library(cowplot)
title <- ggdraw() + draw_label("Posterior joint distribution of Gompertz shape and rate showing true parameter values", fontface='bold')

ggsave(plot = cowplot::plot_grid( title,
                                  cowplot::plot_grid(p1, p2, p3, ncol = 3), ncol = 1, rel_heights=c(0.1, 1) ),

       filename = "pp_joint_gompertz.png", height = 3.8, width = 7, units = "in", dpi = 600)

p <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
p <- cowplot::plot_grid(title, p, ncol = 1, rel_heights=c(0.1, 1))
p
## ----plot-beta-vs-test-----------------------------------------------------
p2 <- ggplot(data.frame(beta = pp_beta01)) +
  geom_density(aes(x = beta)) +
  geom_vline(aes(xintercept = betas01_t[1]), colour = 'red') +
  ggtitle('Posterior distribution of beta\nshowing true value')

cowplot::plot_grid(p1, p2, labels = c("Theta", "beta"))

## ------------------------------------------------------------------------
mean(pp_lambda01 >= lambdas01_t)

## ------------------------------------------------------------------------
mean(pp_gamma01 >= gammas01_t)

## ------------------------------------------------------------------------
mean(pp_lambda01 >= lambdas01_t & pp_gamma01 >= gammas01_t)


### Royston Parmar hazard model

Sys.time()
betas01_t = c(trt = -0.22)
betas02_t = c(trt = -0.4)
betas12_t = c(trt = -0.33)
lambdas01_t = 0.325
lambdas02_t = 0.36
lambdas12_t = 0.39

fn <- function(t, x, betas, ...)
  (-1 + 0.02 * t - 0.03 * t ^ 2 + 0.005 * t ^ 3 + x * betas)
cens = c(4.5, 5.5)

set.seed(9911)
covs <- data.frame(id = 1:2000, trt = stats::rbinom(2000, 1L, 0.5))

sim_rp <- rsimid(
  dist01 = "exponential",
  dist02 = "exponential",
  dist12 = "exponential",
  betas01 = betas01_t,
  betas02 = betas02_t,
  betas12 = betas12_t,
  lambdas01 = lambdas01_t,
  lambdas02 = lambdas02_t,
  lambdas12 = lambdas12_t,
  x = covs,
  cens = cens
)

sim_rp$time_diff <- sim_rp$os_time - sim_rp$df_time

Sys.time()

stanfit3 <- idm_stan(formula01 = Surv(time=df_time,event=df_event)~trt,
                     formula02 = Surv(time=os_time,event=os_event)~trt,
                     formula12 = Surv(time=time_diff,event=os_event)~trt,
                     data = sim_rp,
                     basehaz01 = "ms",
                     basehaz02 = "ms",
                     basehaz12 = "ms",
                     prior01           = rstanarm::normal(),
                     prior_intercept01 = rstanarm::normal(),
                     prior_aux01       = rstanarm::normal(),
                     prior02           = rstanarm::normal(),
                     prior_intercept02 = rstanarm::normal(),
                     prior_aux02       = rstanarm::normal(),
                     prior12           = rstanarm::normal(),
                     prior_intercept12 = rstanarm::normal(),
                     prior_aux12       = rstanarm::normal(),
                     iter = 2000,
                     chains = 4,
                     control = list(adapt_delta = 0.99)
)

stanfit3



print(fit)
rstan::traceplot(fit, 'lp__')

rstan::traceplot(fit, c('beta01','beta02', 'beta12'), ncol = 1)



#### loo

loo1 <- loo(stanfit, cores = 2, k_threshold = 0.7)
loo2 <- loo(stanfit2, cores = 2, k_threshold = 0.7)
loo3 <- loo(stanfit3, cores = 2, k_threshold = 0.7)
compare_models(loo1, loo2, loo3)

#### posterior fit


newdata <- lapply(stanfit2$data, function(x) x[x$id == 13, ])


ps3 <- posterior_fit(stanfit2,
                     newdata = newdata,
                     times = 0,
                     extrapolate = TRUE, condition = FALSE, control = list(edist = 5), type = "cumhaz")

ps3
plot(ps3,
     ids = lapply(seq_along(ps2), function(x) 1:6),
     xlab = list("Time(years)", "Time(years)", "Time since event 1 (years)"),
     labels = c("0 -> 1", "0 -> 2", "1 -> 2"))

#### To conduct the analysis of survival brier score

library(dplyr)
haz <-  as.data.frame(ps3) %>%
  filter(id == 1) %>%
  dplyr::mutate(Haz = median,
                trans = transition) %>%
  select(time, Haz, trans)

library(mstate)
# transition matrix for illness-death model
tmat <- trans.illdeath()

tv <- sort(unique(haz$time))
out <- mssample(Haz=haz,
                trans=tmat,
                clock = "reset",
                M=10000,
                output = "state",
                tvec = tv,
                do.trace=2500)
