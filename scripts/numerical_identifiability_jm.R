library(simms)
# library(rstanarm, lib.loc = "~/R-dev/")
library(rstanarm)
devtools::document()
set.seed(9911)
sim1 = simms::sim_idm_jm(n = 2000)

tmat_mst <- mstate::trans.illdeath(names=c("diagnosis","relapse","death"))
sim1$Mstate = mstate::msprep(time=c(NA,"df_time","os_time"),
                             status=c(NA,"df_event","os_event"),
                             data = sim1$Event,
                             trans=tmat_mst)

sim1$Mstate <- dplyr::left_join(sim1$Mstate,
                                sim1$Event)

formulaLong = as.formula(Yij_1 ~  tij + Z1 + Z2 + tij2 + tij3 + (1 | id))

formulaMs = lapply(1:3, function (x)
  as.formula(Surv(time=time,event=status) ~ Z1 + Z2) )

basehaz = lapply(1:3, function(x)
  "weibull")

priorLong = rstanarm::normal()
priorLong_intercept = rstanarm::normal()
priorLong_aux = rstanarm::cauchy(0, 5)
prior_covariance = rstanarm::lkj()

priorMs_intercept = lapply(1:3, function(x)
  normal() )

priorMs_aux = lapply(1:3, function(x)
  rstanarm::normal() )

priorMs = lapply(1:3, function(x)
  rstanarm::normal() )

priorMs_assoc = lapply(1:3, function(x)
  rstanarm::normal() )

priorEvent_assoc = lapply(1:3, function(x)
  rstanarm::normal() )

dataLong = sim1$Long1
dataLong$tij2 <- dataLong$tij*dataLong$tij  
dataLong$tij3 <- dataLong$tij2*dataLong$tij  
dataMs = sim1$Mstate

time_var = "tij"
time_start = "Tstart"
id_var = "id"

family = gaussian
assoc = "etavalue"
lag_assoc = 0
epsilon = 1E-5
prior_PD = FALSE
qnodes = 15
init = "prefit"
transition_labels = c("01", "02", "12")

options(mc.cores = 4L)
stanfit = msjm_stan(formulaLong = formulaLong,
                    dataLong = dataLong,
                    formulaMs = formulaMs,
                    dataMs = dataMs,
                    time_var = "tij",
                    transition_labels = transition_labels,
                    time_start = "Tstart",
                    id_var = "id",
                    init = "prefit",
                    family = gaussian,
                    assoc = "etavalue",
                    lag_assoc = 0,
                    epsilon = 1E-5,
                    prior_PD = FALSE,
                    priorLong = rstanarm::normal(),
                    priorLong_intercept = rstanarm::normal(),
                    priorLong_aux = rstanarm::cauchy(0, 5),
                    prior_covariance = rstanarm::lkj(),
                    priorMs_intercept = lapply(1:3, function(x)
                      rstanarm::normal() ),
                    priorMs_aux = lapply(1:3, function(x)
                      rstanarm::normal() ),
                    priorMs = lapply(1:3, function(x)
                      rstanarm::normal() ),
                    priorMs_assoc = lapply(1:3, function(x)
                      rstanarm::normal() ),
                    basehaz = lapply(1:3, function(x)
                      "weibull"),
                    iter = 2000, 
                    chains = 4)

saveRDS(nlist(sim1, stanfit), "~/rfactory/mstte-data/jm_stanfit2.RDS")

stanfit = readRDS("~/rfactory/mstte-data/jm_stanfit.RDS")
sim1 = readRDS("~/rfactory/mstte-data/jm_sim1.RDS")
