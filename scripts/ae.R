####### Test on NSCLC data for adverse event modeling
library(mstate)
library(dplyr)
roxygen2::roxygenise(clean = TRUE)

ae = readxl::read_xlsx("~/mstte-data/Data for transmodel.xlsx")
ae$Time2 = ae$Time1 + ae$Time2
ae$Time3 = ae$Time2 + ae$Time3
ae$Time4 = ae$Time4 + ae$Time3

ae$Time6 = ae$Time6 + ae$Time5
ae$Time7 = ae$Time6 + ae$Time7
ae$Time8 = ae$Time8 + ae$Time7
ae$Time10 = ae$Time9 + ae$Time10

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

 ae = dplyr::left_join( ae %>%
  select("ID" ,"Time1", "Time5", "Time9") %>%
  mutate(RowMax = apply(select(.,starts_with("Time")), 1, FUN=max, na.rm=TRUE)) %>%
mutate_at(vars(starts_with("Time")), funs(ifelse(is.na(.), RowMax,.))) %>%
select(-RowMax),
ae %>%
  select(-"Time1", -"Time5", -"Time9") )

 ae = dplyr::left_join( ae %>%
  select("ID" ,"Time1", "Time2", "Time3", "Time4") %>%
  mutate(RowMax = apply(select(.,starts_with("Time")), 1, FUN=max, na.rm=TRUE)) %>%
  mutate_at(vars(starts_with("Time")), funs(ifelse(is.na(.), RowMax,.))) %>%
  select(-RowMax) ,
  ae %>%
    select(-"Time1", -"Time2", -"Time3", -"Time4") )

 ae = dplyr::left_join( ae %>%
                          select("ID" ,"Time5", "Time6", "Time7", "Time8") %>%
                          mutate(RowMax = apply(select(.,starts_with("Time")), 1, FUN=max, na.rm=TRUE)) %>%
                          mutate_at(vars(starts_with("Time")), funs(ifelse(is.na(.), RowMax,.))) %>%
                          select(-RowMax) ,
                        ae %>%
                          select(-"Time5", -"Time6", -"Time7", -"Time8") )

 ae = dplyr::left_join( ae %>%
                          select("ID" ,"Time9", "Time10") %>%
                          mutate(RowMax = apply(select(.,starts_with("Time")), 1, FUN=max, na.rm=TRUE)) %>%
                          mutate_at(vars(starts_with("Time")), funs(ifelse(is.na(.), RowMax,.))) %>%
                          select(-RowMax) ,
                        ae %>%
                          select(-"Time9", -"Time10") )

 status$ID = seq_len(nrow(status))
 ae = left_join(ae, status)

 trans_mat_short =  mstate::transMat( list(
   c(2),
   c(3),
   c(4),
   c(5), c() ), c("START", "LRI1", "RLRI1", "LRI2", "RLRI2"))

 ae_mstate1 <- mstate::msprep(
   time=c(NA,"Time5",  "Time6", "Time7", "Time8"),
   status=c(NA,"status5",  "status6", "status7", "status8"),
   data= cbind(ae, status),
   trans=trans_mat_short )

 trans_mat_short =  mstate::transMat( list(
   c(2),
   c(3),
   c(4),
   c(5), c() ), c("START", "DYS1", "RDYS1", "DYS2", "RDYS2"))

 ae_mstate2 <- mstate::msprep(
   time=c(NA,"Time1",  "Time2", "Time3", "Time4"),
   status=c(NA,"status1",  "status2", "status3", "status4"),
   data= ae,
   trans=trans_mat_short )
 ae_mstate2$trans = ae_mstate2$trans + max(ae_mstate1$trans)

 trans_mat_short =  mstate::transMat( list(
   c(2),
   c(3),
   c() ), c("START", "ILDLRI", "RILDLRI"))

 ae_mstate3 <- mstate::msprep(
   time=c(NA,"Time9",  "Time10"),
   status=c(NA,"status9",  "status10"),
   data= ae,
   trans=trans_mat_short )
 ae_mstate3$trans = ae_mstate3$trans + max(ae_mstate2$trans)


 ae_mstate = rbind(ae_mstate1, ae_mstate2, ae_mstate3)

 ae_mstate = ae_mstate[order(ae_mstate$id, ae_mstate$time), ]

 colnames(ae) <- tolower(colnames(ae))
 ae_mstate <- left_join(ae_mstate, ae[ ,c("id", "age", "sex")])

formula = lapply(1:max(ae_mstate$trans), function (x)
  as.formula(Surv(time=time,event=status)~age+sex) )
# prior = lapply(1:3, function(x)
#   rstanarm::student_t() )

prior_intercept = lapply(1:max(ae_mstate$trans), function(x)
  rstanarm::normal() )

prior_aux = lapply(1:max(ae_mstate$trans), function(x)
  rstanarm::normal() )

basehaz = lapply(1:max(ae_mstate$trans), function(x)
  "ms")

options(mc.cores = 2)
stanfit <- mstte_stan(formula = formula,
                      data = ae_mstate %>%
                        filter(time > 0),
                      # transition_labels = c("DP", "DX", "DPDX"),
                      basehaz = basehaz,
                      prior           = lapply(1:max(ae_mstate$trans), function(x)
                        rstanarm::normal() ),
                      prior_intercept = prior_intercept,
                      prior_aux       = prior_aux,
                      iter = 1000,
                      chains = 4,
                      control = list(adapt_delta = 0.95)
)


