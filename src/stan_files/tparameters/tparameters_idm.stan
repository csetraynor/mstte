
  vector[e_K01] e_beta01;               // log hazard ratios
  vector[e_K02] e_beta02;               // log hazard ratios
  vector[e_K12] e_beta12;               // log hazard ratios

  vector[a_K01] a_beta01;               // assoc params
  vector[a_K02] a_beta02;               // assoc params
  vector[a_K12] a_beta12;               // assoc params

  vector[basehaz_nvars01] e_aux01;      // basehaz params
  vector[basehaz_nvars02] e_aux02;      // basehaz params
  vector[basehaz_nvars12] e_aux12;      // basehaz params
 //---- Parameters for event submodel
  e_beta01 = make_beta(e_z_beta01, e_prior_dist01, e_prior_mean01,
                     e_prior_scale01, e_prior_df01, e_global_prior_scale01,
                     e_global01, e_local01, e_ool01, e_mix01, rep_array(1.0, 0), 0,
                     e_slab_scale01, e_caux01);
e_beta02 = make_beta(e_z_beta02, e_prior_dist02, e_prior_mean02,
                     e_prior_scale02, e_prior_df02, e_global_prior_scale02,
                     e_global02, e_local02, e_ool02, e_mix02, rep_array(1.0, 0), 0,
                     e_slab_scale02, e_caux02);

e_beta12 = make_beta(e_z_beta12, e_prior_dist12, e_prior_mean12,
                     e_prior_scale12, e_prior_df12, e_global_prior_scale12,
                     e_global12, e_local12, e_ool12, e_mix12, rep_array(1.0, 0), 0,
                     e_slab_scale12, e_caux12);


  a_beta01 = make_beta(a_z_beta01, a_prior_dist01, a_prior_mean01,
                     a_prior_scale01, a_prior_df01, a_global_prior_scale01,
                     a_global01, a_local01, a_ool01, a_mix01, rep_array(1.0, 0), 0,
                     a_slab_scale01, a_caux01);
  a_beta02 = make_beta(a_z_beta02, a_prior_dist02, a_prior_mean02,
                     a_prior_scale02, a_prior_df02, a_global_prior_scale02,
                     a_global02, a_local02, a_ool02, a_mix02, rep_array(1.0, 0), 0,
                     a_slab_scale02, a_caux02);
  a_beta12 = make_beta(a_z_beta12, a_prior_dist12, a_prior_mean12,
                     a_prior_scale12, a_prior_df12, a_global_prior_scale12,
                     a_global12, a_local12, a_ool12, a_mix12, rep_array(1.0, 0), 0,
                     a_slab_scale12, a_caux12);

  e_aux01  = make_basehaz_coef(e_aux_unscaled01, e_prior_dist_for_aux01,
                             e_prior_mean_for_aux01, e_prior_scale_for_aux01);
  e_aux02  = make_basehaz_coef(e_aux_unscaled02, e_prior_dist_for_aux02,
          e_prior_mean_for_aux02, e_prior_scale_for_aux02);
  e_aux12  = make_basehaz_coef(e_aux_unscaled12, e_prior_dist_for_aux12,
                             e_prior_mean_for_aux12, e_prior_scale_for_aux12);
