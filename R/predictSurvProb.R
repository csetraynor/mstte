predictSurvProb <- function(object,newdata,times,...){
  ## baselineHazard.coxph(object,times)
  ## require(survival)
  ## new feature of the survival package requires that the
  ## original data are included
  ## survival.survfit.coxph <- getFromNamespace("survfit.coxph",ns="survival")
  ## survival.summary.survfit <- getFromNamespace("summary.survfit",ns="survival")
  ## b <- function(x){browser()}
  ## b()
  survfit.object <- survival::survfit(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
  if (is.null(attr(object$terms,"specials")$strata)){
    ## case without strata
    inflated.pred <- summary(survfit.object,times=times)$surv
    if (is.null(inflated.pred)){
      ## can happen when all times beyond maxtime
      p=matrix(NA,ncol=length(times),nrow=NROW(newdata))
    } else{
      p <- t(inflated.pred)
      if ((beyond <- (length(times)-NCOL(p)))>0)
        p <- cbind(p,matrix(NA,nrow=NROW(newdata),ncol=beyond))
    }
  }else{
    ## case with strata
    inflated.pred <- summary(survfit.object,times=times)
    plist <- split(inflated.pred$surv,inflated.pred$strata)
    p <- do.call("rbind",lapply(plist,function(x){
      beyond <- length(times)-length(x)
      c(x,rep(NA,beyond))
    }))
    ## p <- matrix(inflated.pred,ncol=length(times),byrow=TRUE)
  }
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}
