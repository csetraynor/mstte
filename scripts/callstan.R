library(rstanarm)
library(rstan)
devtools::document()
options(mc.cores = 2)
rstan_options(auto_write = TRUE)

library(mstate)
data("ebmt1")

tmat.mst <- trans.illdeath(names=c("transplant","relapse","death"))
tmat.mst <- matrix(NA,3,3)
tmat.mst[1,2:3] <- 1:2
tmat.mst[2,3] <- 3
dimnames(tmat.mst) <- list(c("transplant","relapse","death"),
                           c("transplant","relapse","death"))

ebmt.mstate <- msprep(time=c(NA,"rel","srv"), status=c(NA,"relstat","srvstat"), data=ebmt1, trans=tmat.mst , keep = c("age", "score"))

survreg(Surv(time=time,event=status)~strata(trans) * age + strata(trans) * score, data = ebmt.mstate, dist = "wei")

formula = lapply(1:3, function (x)
  as.formula(Surv(time=time,event=status)~score) )
formula[[1]] =  as.formula(Surv(time=time,event=status)~score)
formula[[2]] =  as.formula(Surv(time=time,event=status)~age)


library(dplyr)
library(survival)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
stanfit <- idm_stan(formula01 = Surv(time=rel,event=relstat)~score,
                    formula02 = Surv(time=srv,event=srvstat)~age,
                    formula12 = Surv(time=time_diff,event=srvstat)~score+age,
                    data = ebmt1 %>%
                      mutate(time_diff = srv - rel),
                    basehaz01 = "weibull",
                    basehaz02 = "exp",
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

# prior = lapply(1:3, function(x)
#   rstanarm::student_t() )

prior_intercept = lapply(1:3, function(x)
  rstanarm::normal() )

prior_aux = lapply(1:3, function(x)
  rstanarm::normal() )

basehaz = lapply(1:3, function(x)
  "ms")
basehaz[[1]] = "weibull"
basehaz[[2]] = "exp"


stanfit <- mstte_stan(formula = formula,
                     data = ebmt.mstate,
                     transition_labels = c("DP", "DX", "DPDX"),
                     basehaz = basehaz,
                     prior           = lapply(1:3, function(x)
                       rstanarm::normal() ),
                     prior_intercept = prior_intercept,
                     prior_aux       = prior_aux,
                     iter = 1000,
                     chains = 4,
                     control = list(adapt_delta = 0.95)
)

coxph(Surv(time=time,event=status)~age+score+strata(trans), data = ebmt.mstate)


stanfit


library(mstate)
data("aidssi2")

trans_mat = transMat(
  x = list(c(2, 3), c(4, 5), c(4), c(5), c()),
  names = c("HIV", "AIDS", "SI", "AIDS/SI", "death") )

aidssi2 <- within(aidssi2, {
  siaids.time = pmax(aids.time,si.time)
  siaids.stat = ifelse(aids.stat & si.stat, 1, 0)
}
       )

aidss_mst <- msprep(
  time = c(NA, "aids.time", "si.time", "siaids.time", "death.time"),
  status = c(NA, "aids.stat", "si.stat", "siaids.stat", "death.stat") ,
  data = aidssi2, trans = trans_mat,
  start = list(time = aidssi2$entry.time, state = rep(1, nrow(aidssi2))),
  id = "patnr"
)

aidss_mst <- dplyr::left_join(aidss_mst, aidssi2[ ,c("patnr", "age.inf", "ccr5")])


formula = lapply(1:6, function (x)
  as.formula(Surv(time=time,event=status)~age.inf+ccr5) )
prior_intercept = lapply(1:6, function(x)
  rstanarm::normal() )

prior_aux = lapply(1:6, function(x)
  rstanarm::normal() )

basehaz = lapply(1:6, function(x)
  "ms")
basehaz[3:5] = lapply(1:3, function(x)
  "weibull")
basehaz[[6]] = "exp"

options(mc.cores = 2)
stanfit <- mstte_stan(formula = formula,
                      data = aidss_mst,
                     # transition_labels = c("DP", "DX", "DPDX"),
                      basehaz = basehaz,
                      prior           = lapply(1:6, function(x)
                        rstanarm::normal() ),
                      prior_intercept = prior_intercept,
                      prior_aux       = prior_aux,
                      iter = 1000,
                      chains = 2,
                      control = list(adapt_delta = 0.95)
)

stanfit
summary(stanfit)
loo1 <- loo(stanfit, cores = 2, k_threshold = 0.7)


