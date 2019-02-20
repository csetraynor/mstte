  int<lower=0> e_hs01 = get_nvars_for_hs(e_prior_dist01);
  int<lower=0> a_hs01 = get_nvars_for_hs(a_prior_dist01);
  int<lower=0> e_hs02 = get_nvars_for_hs(e_prior_dist02);
  int<lower=0> a_hs02 = get_nvars_for_hs(a_prior_dist02);
  int<lower=0> e_hs12 = get_nvars_for_hs(e_prior_dist12);
  int<lower=0> a_hs12 = get_nvars_for_hs(a_prior_dist12);

  vector[len_epts01] log_epts01  = log(epts01); // log of event times
  vector[len_epts02] log_epts02  = log(epts02); // log of event times
  vector[len_epts12] log_epts12  = log(epts12); // log of event times

  vector[len_qpts01] log_qpts01  = log(qpts01); // log of quadrature points
  vector[len_qpts02] log_qpts02  = log(qpts02); // log of quadrature points
  vector[len_qpts12] log_qpts12  = log(qpts12); // log of quadrature points

  vector[len_ipts01] log_ipts01  = log(ipts01); // log of qpts for interval censoring
  vector[len_ipts02] log_ipts02  = log(ipts02); // log of qpts for interval censoring
  vector[len_ipts12] log_ipts12  = log(ipts12); // log of qpts for interval censoring

  real sum_epts01     = sum(epts01);     // sum of event times
  real sum_epts02     = sum(epts02);     // sum of event times
  real sum_epts12     = sum(epts12);     // sum of event times

  real sum_log_epts01 = sum(log_epts01); // sum of log event times
  real sum_log_epts02 = sum(log_epts02); // sum of log event times
  real sum_log_epts12 = sum(log_epts12); // sum of log event times
