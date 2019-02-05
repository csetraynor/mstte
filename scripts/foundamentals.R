library(mstate)
data("ebmt1")
data("aidssi")
head(ebmt1)
nrow(ebmt1)
summary(ebmt1)

ebmt1 <- within(ebmt1, {
  time <- pmin(srv,rel)/365.25
  stat <- ifelse(rel<srv,relstat,srvstat*2)
  type <- factor(stat,
                 labels = c("Event-Free","Relapse","Death"))
})


fit <- survfit(Surv(time, stat == 1)~1, data = ebmt1)
summary(fit)
plot(fit, fun = "cumhaz")

fitcox <- coxph(Surv(time, stat == 1)~score, data = ebmt1)
summary(fitcox)
cox.zph(fitcox)

lr <- survdiff( Surv(time, stat == 1)~score, data = ebmt1)
plot

### Competing risks non parametric

csiSurv <- survfit(Surv(time = time, event = status, type = "mstate")~1,
                   data = aidssi)
summary(csiSurv, times=seq(2,10,by=2))
plot(csiSurv)

cuminc.1 <- survfit(Surv(time = time, event = stat, type = "mstate")~1,
                   data = ebmt1)

summ.1 <- summary(cuminc.1, times=c(1,5))
summ.1$lower
summ.1

library(etm)
data(abortion)

cif.ab <- etmCIF(survival::Surv(entry, exit, cause != 0) ~ group, abortion,
                 etype = cause, failcode = 3)

cif.ab

plot(cif.ab, ci.type = "bars", pos.ci = 24,
     col = c(1, 2), lty = 1, curvlab = c("Control", "Exposed"))


## Multi-state model
library(mstate)
data("aidssi2")
head(aidssi2)

trans.mat <- transMat(
  x = list(c(2,3), c(4,5), c(4), c(5), c()),
  names = c("HIV", "AIDS", "SI", "AIDS/SI", "death")
)


aidssi2 <- within(aidssi2, {
  siaids.time <- pmax(si.time, aids.time)
  siaids.stat <- ifelse(aids.stat & si.stat, 1, 0)
})


trans.mat
aidssi.mst <- msprep(
  time = c(NA, "aids.time", "si.time", "siaids.time", "death.time"),
  status = c(NA, "aids.stat", "si.stat", "siaids.stat", "death.stat"),
  data = aidssi2, trans = trans.mat,
  start = list(time = aidssi2$entry.time, state = rep(1, nrow(aidssi2)) )
)



data("ebmt1")

tmat.mst <- trans.illdeath(names=c("transplant","relapse","death"))
tmat.mst <- matrix(NA,3,3)
tmat.mst[1,2:3] <- 1:2
tmat.mst[2,3] <- 3
dimnames(tmat.mst) <- list(c("transplant","relapse","death"),
                           c("transplant","relapse","death"))

ebmt.mstate <- msprep(time=c(NA,"rel","srv"), status=c(NA,"relstat","srvstat"), data=ebmt1, trans=tmat.mst , keep = c("age", "score"))
ebmt1[1:2, ]
head(ebmt1)
ebmt.mstate[1:10, ]
subset(ebmt.mstate, id %in% c(1,2))

scl <- sas7bdat::read.sas7bdat("/media/mtr/A5C2-009E/SCL/c9732_demographic.sas7bdat")

scl <- within(scl, {
  PD_STATUS <- ifelse(is.nan(PD_TIME), 0, 1)
  OS_STATUS <- ifelse(STATUS == 2, 1, 0)
  PD_TIME = PFS_TIME
  TIMEDIFF = OS_TIME - PD_TIME
  TRTARM = as.factor(TRT_ARM)
})

library(mstate)
tmat.mst <- trans.illdeath(names=c("healthy","relapse","death"))
# tmat.mst <- matrix(NA,3,3)
# tmat.mst[1,2:3] <- 1:2
# tmat.mst[2,3] <- 3
# dimnames(tmat.mst) <- list(c("transplant","relapse","death"),
#                            c("transplant","relapse","death"))

scl = scl[scl$PD_TIME > 0, ]

scl.ms <- msprep(time=c(NA,"PD_TIME","OS_TIME"), status=c(NA,"PD_STATUS","OS_STATUS"), data=scl, trans=tmat.mst , keep = "TRTARM")

scl.ms$time[scl.ms$time == 0] <- 1e-10 # add small lag time

formula = lapply(1:3, function (x)
  as.formula(Surv(time=time,event=status)~age+score) )
formula[[1]] =  as.formula(Surv(time=time,event=status)~score)
formula[[2]] =  as.formula(Surv(time=time,event=status)~age)


prior = lapply(1:3, function(x)
  rstanarm::normal() )

prior_intercept = lapply(1:3, function(x)
  rstanarm::normal() )

prior_aux = lapply(1:3, function(x)
  rstanarm::normal() )

basehaz = lapply(1:3, function(x)
  "ms")
basehaz[[1]] = "weibull"
basehaz[[2]] = "exp"

