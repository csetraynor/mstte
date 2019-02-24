library(simms)
library(rstanarm, lib.loc = "~/R-dev/")
#library(rstanarm)
devtools::document()
set.seed(9911)
sim1 = simms::sim_idm_jm(n = 100,
                         seed = 9911,
                         fixed_trajectory = "linear",
                         random_trajectory = "none",
                         b_sd = 1.5,
                         assoc = "etavalue",
                         basehaz = "weibull")

tmat_mst <- mstate::trans.illdeath(names=c("diagnosis","relapse","death"))
sim1$Mstate = mstate::msprep(time=c(NA,"df_time","os_time"),
                             status=c(NA,"df_event","os_event"),
                             data = sim1$Event,
                             trans=tmat_mst)

sim1$Mstate <- dplyr::left_join(sim1$Mstate,
                                sim1$Event)

formulaLong = as.formula(Yij_1 ~  tij + Z1 + Z2 + (1 | id))

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
# dataLong$tij2 <- dataLong$tij*dataLong$tij  
# dataLong$tij3 <- dataLong$tij2*dataLong$tij  
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
                    iter = 1, 
                    chains = 1)

saveRDS(nlist(sim1, stanfit), "~/rfactory/mstte-data/jm_stanfit16.RDS")

## end of simulation

simfit = readRDS("~/rfactory/mstte-data/jmMS_ni2000.RDS")
devtools::document()
sumfit = as.data.frame(summary(simfit$stanfit))
actual = unlist( attr(simfit$sim1,"params") )

sumb = sumfit[grepl("b\\[Long1", rownames(sumfit)), ]
sumfit = sumfit[!grepl("b\\[Long1|posterior", rownames(sumfit)), ]

actual <- c(
  "(Intercept)_1"  = log(lambdas01_t),
  "weibull-shape_1" = gammas01_t,
  
  "(Intercept)_2" = log(lambdas02_t),
  "weibull-shape_2" = gammas02_t,
  
  "(Intercept)_3" = log(lambdas12_t),
  "weibull-shape_3" = gammas12_t)


actual = actual[ c(1,5,2,3,grep("betaEvent_intercept", names(actual) ), 7,8,12,13,17,18,9,14,19,4,grep("betaEvent_aux", names(actual) ),  21, 22) ]

numtest = cbind(actual, sumfit)

posterior <- as.array(simfit$stanfit$stanfit)


