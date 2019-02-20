functions {
#include /functions/common_functions.stan
#include /functions/bernoulli_likelihoods.stan
#include /functions/binomial_likelihoods.stan
#include /functions/continuous_likelihoods.stan
#include /functions/count_likelihoods.stan
#include /functions/mvmer_functions.stan
#include /functions/hazard_functions.stan
#include /functions/jm_functions.stan
}
data {
  // declares:
  //   M
  //   resp_type
  //   intercept_type
  //   has_{aux,weights}
  //   y{Nobs,Neta,K}
  //   b{N1,N2}
  //   b{K1,K2}
  //   b{K1_len,K2_len}
  //   b{K1_idx,bK2_idx}
  //   t, p, l, q
  //   len_theta_L
#include /data/dimensions_mvmer.stan

  // declares:
  //   yInt{1,2,3}
  //   yReal{1,2,3}
  //   yX{1,2,3}
  //   yXbar{1,2,3}
  //   y{1,2,3}_Z{1,2}
  //   y{1,2,3}_Z{1,2}_id
  //   y_prior_dist{_for_intercept,_for_aux,_for_cov}
  //   family
  //   link
  //   prior_PD
#include /data/data_mvmer.stan

  // declares:
  //   e_prior_dist{_for_intercept,_for_aux}
  //   e_{K,x,basis,has_intercept}
  //   qwts_{event,lcens,rcens,icenl,icenu,delay}
  //   q{event,lcens,rcens,icenl,icenu,delay}
  //   qnodes
  //   cpts
  //   len_cpts
  //   idx_cpts
  //   basehaz_type
  //   norm_const
#include /data/data_event01.stan
#include /data/data_event02.stan
#include /data/data_event12.stan

  // declares:
  //   a_{K,xbar}
  //   a_prior_dist, assoc, assoc_uses, has_assoc
  //   {sum_}size_which_b, which_b_zindex
  //   {sum_}size_which_coef, which_coef_{zindex,xindex}
  //   {sum_,sum_size_}which_interactions
  //   y_nrow{_auc}_cpts
  //   auc_{qnodes,qwts}
  //   a_K_data
  //   has_grp
  //   grp_assoc
  //   idx_grp
  //   idx_cpts
  //   y{1,2,3}_x_{eta,eps,auc}_cpts
  //   y{1,2,3}_z{1,2}_{eta,eps,auc}_cpts
  //   y{1,2,3}_z{1,2}_id_{eta,eps,auc}_cpts
#include /data/data_assoc01.stan
#include /data/data_assoc02.stan
#include /data/data_assoc12.stan

  // declares:
  //   e_prior_{mean,scale,df}{_for_intercept,for_aux},
  //   e_global_prior_{scale,df}
#include /data/hyperparameters_mvmer.stan
#include /data/hyperparameters_idm.stan

  // declares:
  //   a_prior_{mean,scale,df},
  //   a_global_prior_{scale,df}
#include /data/hyperparameters_assoc_idm.stan
}

transformed data {
  //declares
  // e_hs{01,02,12}
  // a_hs{01,02,12}
  // log_epts{01,02,12}
  // log_ipts{01,02,12}
#include /tdata/tdata_idm.stan

  // declares:
  //   yHs{1,2,3}
  //   len_{z_T,var_group,rho}
  //   pos
  //   delta
  //   bCov{1,2}_idx
  //   {sqrt,log,sum_log}_y{1,2,3}
#include /tdata/tdata_mvmer.stan

}

parameters {
  // declares:
  //   yGamma{1,2,3}
  //   z_yBeta{1,2,3}
  //   z_b
  //   z_T
  //   rho
  //   zeta
  //   tau
  //   bSd{1,2}
  //   z_bMat{1,2}
  //   bCholesky{1,2}
  //   yAux{1,2,3}_unscaled
  //   yGlobal{1,2,3}
  //   yLocal{1,2,3},
  //   yOol{1,2,3}
  //   yMix{1,2,3}
#include /parameters/parameters_mvmer.stan

  // declares:
  //   e_{gamma,z_beta,aux_unscaled,global,local,mix,ool}
#include /parameters/parameters_idm.stan
  // declares:
  //   a_{z_beta,global,local,mix,ool}
#include /parameters/parameters_assoc_idm.stan
}

transformed parameters {
  // declares and defines:
  //   yBeta{1,2,3}
  //   yAux{1,2,3}
  //   yAuxMaximum,
  //   theta_L
  //   bMat{1,2}
#include /tparameters/tparameters_mvmer.stan

  // declares and defines:
  //  e_beta{01,02,12}
  //  a_beta{01,02,12}
  //  e_aux{01,02,12}
#include /tparameters/tparameters_idm.stan
}

model {

  //---- Log likelihoods for longitudinal submodels
#include /model/mvmer_lp.stan


  //---- Log likelihood for event submodel (GK quadrature)
  {
    vector[len_cpts01] e_eta01;
    vector[len_cpts02] e_eta02;
    vector[len_cpts12] e_eta12;

    // Event submodel: linear predictor at event and quad times
    if (e_K01 > 0) e_eta01 = e_x01 * e_beta01;
    else         e_eta01 = rep_vector(0.0, len_cpts01);

    if (e_K01 > 0) e_eta02 = e_x02 * e_beta02;
    else         e_eta02 = rep_vector(0.0, len_cpts02);

    if (e_K01 > 0) e_eta12 = e_x12 * e_beta12;
    else         e_eta12 = rep_vector(0.0, len_cpts12);


    if (e_has_intercept01 == 1) e_eta01 = e_eta01 + e_gamma01[1];
    if (e_has_intercept02 == 1) e_eta02 = e_eta02 + e_gamma02[1];
    if (e_has_intercept12 == 1) e_eta12 = e_eta12 + e_gamma12[1];

    if (assoc01 == 1) {
      // declares:
      //   y_eta_q{_eps,_lag,_auc}
      //   y_eta_qwide{_eps,_lag,_auc}
      //   y_q_wide{_eps,_lag,_auc}
      //   mark{2,3}
#include /model/assoc_evaluate01.stan
    }

    if (assoc02 == 1) {
      // declares:
      //   y_eta_q{_eps,_lag,_auc}
      //   y_eta_qwide{_eps,_lag,_auc}
      //   y_q_wide{_eps,_lag,_auc}
      //   mark{2,3}
#include /model/assoc_evaluate02.stan
    }

    if (assoc12 == 1) {
      // declares:
      //   y_eta_q{_eps,_lag,_auc}
      //   y_eta_qwide{_eps,_lag,_auc}
      //   y_q_wide{_eps,_lag,_auc}
      //   mark{2,3}
#include /model/assoc_evaluate12.stan
    }

      {
    // declares (and increments target with event log-lik):
    //   log_basehaz,
    //   log_{haz_q,haz_etimes,surv_etimes,event}
#include /model/event_lp01.stan
#include /model/event_lp02.stan
#include /model/event_lp12.stan
    }

  }

  //---- Log priors
  // increments target with mvmer priors
#include /model/priors_mvmer.stan
  // increments target with idm priors
#include /model/priors_idm.stan
}

generated quantities {
  // declares and defines:
  //   mean_PPD
  //   yAlpha{1,2,3}
  //   b{1,2}
  //   bCov{1,2}
#include /gqs/gen_quantities_mvmer.stan

  // declares and defines:
  // e_alpha{01,02,12}

  // norm_const is a constant shift in log baseline hazard
#include /gqs/gen_quantities_idm.stan
}

