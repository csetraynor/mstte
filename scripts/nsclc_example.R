library(dplyr)
require(sas7bdat)

tte = sas7bdat::read.sas7bdat("~/rfactory/mstte-data/NSCLCstageIV/jfcc_deid_ctrl_trt/adtte.sas7bdat")
dplyr::glimpse(tte)
unique(tte$PARAM)

pfs = tte %>%
  filter(PARAMCD == "PFS") %>%
  select("USUBJID", pfs = "EVNTDESC")
pfs1 = tte[tte$PARAMCD == "PFSS1", c("USUBJID", "EVNTDESC")]
pfs2 = tte[tte$PARAMCD == "PFSS2", c("USUBJID", "EVNTDESC")]
pfs3 = tte[tte$PARAMCD == "PFSS3", c("USUBJID", "EVNTDESC")]
pfs4 = tte[tte$PARAMCD == "PFSS4", c("USUBJID", "EVNTDESC")]

PFS = Reduce(function(x,y) merge(x = x, y = y, by = "USUBJID"), 
       list(tte %>%
              filter(PARAMCD == "PFS") %>%
              select("USUBJID", pfs = "EVNTDESC"),
            tte %>%
              filter(PARAMCD == "PFSS1") %>%
              select("USUBJID", pfss1 = "EVNTDESC"),
            tte %>%
              filter(PARAMCD == "PFSS2") %>%
              select("USUBJID", pfss2 = "EVNTDESC"),
            tte %>%
              filter(PARAMCD == "PFSS3") %>%
              select("USUBJID", pfss3 = "EVNTDESC"),
            tte %>%
              filter(PARAMCD == "PFSS4") %>%
              select("USUBJID", pfss4 = "EVNTDESC")) )

pfs =  tte %>%
  filter(PARAMCD == "PFSS2") %>%
  select("USUBJID", event_desc_pd = "EVNTDESC", dp_time = "ADY") %>%
  mutate(event_desc_pd = as.character(event_desc_pd),
         dp_status = event_desc_pd == "Disease Progression")

os = tte %>%
  filter(PARAMCD == "OS") %>%
  select("USUBJID", event_desc_os = "EVNTDESC", os_time = "ADY") %>%
  mutate(event_desc_os = as.character(event_desc_os),
         os_status = event_desc_os == "Death")

nscl_surv = left_join(os, pfs, by = "USUBJID")


### Add covariates to survival submodel(s)
sl = sas7bdat::read.sas7bdat("~/rfactory/mstte-data/NSCLCstageIV/jfcc_deid_ctrl_trt/adsl.sas7bdat")
table(sl$EGFKTIFL)

sl = sl %>%
  select(USUBJID, 
         BMIBL, AGE, ECOGBL, SEX,
         EGFKTIFL,  
         COUNGR2,SMKCAT1N )

nscl_surv <- left_join(nscl_surv, sl, by =  "USUBJID")
sapply(nscl_surv, function(x) sum(is.na(x)) )

nscl_surv$BMIBL[is.na(nscl_surv$BMIBL)] <- median(nscl_surv$BMIBL, na.rm = TRUE)

## Longitudinal submodels

lab = sas7bdat::read.sas7bdat("~/rfactory/mstte-data/NSCLCstageIV/jfcc_deid_ctrl_trt/adlb.sas7bdat")
dplyr::glimpse(lab)
unique(lab$LBTEST)

# LDH biomarker

ldh = lab[lab$LBTESTCD == "LDH", ]
dplyr::glimpse(ldh)

# Parameter value in U/L
nscl_ldh = ldh %>%
  select(USUBJID,
         time = ADY, 
         ldh = AVAL, 
         PARAMCD,
         PARAMUNT, PARCAT1) %>%
  filter(!is.na(ldh))

library(data.table)

nsclc_ldh_baseline <- data.table(nscl_ldh, key = "USUBJID")
nsclc_ldh_baseline <- nsclc_ldh_baseline[, head(.SD, 1), by = key(nsclc_ldh_baseline)]

nscl_surv <- left_join(nscl_surv,
                       nsclc_ldh_baseline %>%
                         mutate(ldhBL = ldh) %>%
                         select(ldhBL, USUBJID))

nscl_ldh = nscl_ldh %>%
  filter(time >= 0)



nscl_ldh$USUBJID <- as.character(nscl_ldh$USUBJID)
nscl_surv$USUBJID <- as.character(nscl_surv$USUBJID)
nscl_surv <- nscl_surv[ nscl_surv$USUBJID %in% nscl_ldh$USUBJID, ]

nscl_surv$id <- seq_along(nscl_surv$USUBJID )

nscl_ldh = left_join(nscl_ldh, nscl_surv, by = "USUBJID")
sapply(nscl_ldh, function(x) sum(is.na(x)) )


## Monocyte counts
mono = lab[lab$LBTEST == "Monocytes", ]




saveRDS(nlist(nscl_ldh, nscl_surv), "~/rfactory/mstte-data/nscl_data.RDS" )
