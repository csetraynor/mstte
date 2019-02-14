library(simms)

set.seed(9911)
sim1 = simms::sim_idm_jm(n = 2000)

tmat_mst <- mstate::trans.illdeath(names=c("diagnosis","relapse","death"))
sim1$Mstate = mstate::msprep(time=c(NA,"df_time","os_time"),
                             status=c(NA,"df_event","os_event"),
                             data = sim1$Event,
                             trans=tmat_mst)

sim1$Mstate <- dplyr::left_join(sim1$Mstate,
                                sim1$Event)

formulaLong = as.formula(Yij_1 ~  tij + (1 | id))

formulaMs = lapply(1:3, function (x)
  as.formula(Surv(time=time,event=status) ~ Z1 + Z2) )

basehaz = lapply(1:3, function(x)
  "weibull")

prior_intercept = lapply(1:3, function(x)
  rstanarm::normal() )

prior_aux = lapply(1:3, function(x)
  rstanarm::normal() )
