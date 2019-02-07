
library(rstanarm)
library(rstan)
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
