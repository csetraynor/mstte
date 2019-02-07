####### Test on NSCLC data for adverse event modeling
library(mstate)
library(dplyr)
roxygen2::roxygenise(clean = TRUE)

# ae = readxl::read_xlsx("~/rfactory/mstte-data/Data for transmodel.xlsx")
# #
# ae$Time2 = ae$Time1 + ae$Time2
# ae$Time3 = ae$Time2 + ae$Time3
# ae$Time4 = ae$Time4 + ae$Time3
#
# ae$Time6 = ae$Time6 + ae$Time5
# ae$Time7 = ae$Time6 + ae$Time7
# ae$Time8 = ae$Time8 + ae$Time7
# ae$Time10 = ae$Time9 + ae$Time10
# saveRDS(ae, "~/rfactory/mstte-data/ae_datatest.RDS")
ae = readRDS("~/rfactory/mstte-data/ae_datatest.RDS")

nt = 10

ae_time = ae[ ,grepl("Time", colnames(ae))]
status = lapply(seq_len(nt), function(i){
  tmp = as.integer(!is.na(ae_time[i]))
  out = assign(
    paste("status", i, sep = ""), tmp
  )
  out
}
  )
status = do.call(cbind.data.frame, status)
colnames(status) = paste("status", seq_len(ncol(status)), sep = "" )

ae = ae %>% mutate(RowMax = apply(select(.,starts_with("Time")), 1, FUN=max, na.rm=TRUE)) %>%
mutate_at(vars(starts_with("Time")), funs(ifelse(is.na(.), RowMax,.))) %>%
select(-RowMax)

ae = cbind(ae, status)

trans_mat_mstate = mstate::transMat( list(
  c(2,6,10),
  c(3),
  c(4),
  c(5),
  c(),
  c(7),
  c(8),
  c(9),
  c(),
  c(11),
  c()
) ,
c("START", "DYS1", "RDYS1", "DYS2", "RDYS2", "LRI1", "RLRI1", "LRI2", "RLRI2", "ILDLRI", "RILDLRI"))

trans_mat_etm = trans_mat_mstate
trans_mat_etm[is.na(trans_mat_etm)] <- FALSE
trans_mat_etm[trans_mat_etm != 0 ] <- TRUE

library(etm)


head(ae)

ae_mstate <- mstate::msprep(
  time=c(NA,"Time1","Time5", "Time9", "Time2", "Time3", "Time4", "Time6", "Time7", "Time8", "Time10"),
  status=c(NA,"status1","status5", "status9", "status2", "status3", "status4", "status6", "status7", "status8", "status10"),
  data=ae,
  trans=trans_mat_mstate )

ae_etm = etm::etmprep( time=c(NA,"Time1","Time5", "Time9", "Time2", "Time3", "Time4", "Time6", "Time7", "Time8", "Time10"),
                       status=c(NA,"status1","status5", "status9", "status2", "status3", "status4", "status6", "status7", "status8", "status10"),
                       data=ae, cens.name = "cens",
                       tra= trans_mat_etm
)
ae_etm <- ae_etm %>%
  mutate(cens = ifelse(to == "cens", 0, 1),
         time = exit - entry) %>%
  filter(to != "cens")

ae_etm$trans <-  as.numeric( factor(paste(ae_etm$from, ae_etm$to)) )

ae_etm <- ae_etm[ae_etm$time > 0 , ]

colnames(ae) <- tolower(colnames(ae))
ae_etm <- left_join(ae_etm, ae[ ,c("id", "age", "sex")])

formula = lapply(1:max(ae_etm$trans), function (x)
  as.formula(Surv(time=time,event=cens)~age+sex) )
# prior = lapply(1:3, function(x)
#   rstanarm::student_t() )

prior_intercept = lapply(1:max(ae_etm$trans), function(x)
  rstanarm::normal() )

prior_aux = lapply(1:max(ae_etm$trans), function(x)
  rstanarm::normal() )

basehaz = lapply(1:max(ae_etm$trans), function(x)
  "weibull")
basehaz[[1]] = "weibull"
basehaz[[2]] = "exp"


stanfit <- mstte_stan(formula = formula,
                      data = ae_etm,
                      # transition_labels = c("DP", "DX", "DPDX"),
                      basehaz = basehaz,
                      prior           = lapply(1:max(ae_etm$trans), function(x)
                        rstanarm::normal() ),
                      prior_intercept = prior_intercept,
                      prior_aux       = prior_aux,
                      iter = 1000,
                      chains = 4,
                      control = list(adapt_delta = 0.95)
)


