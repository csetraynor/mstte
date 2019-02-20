beta_lp(e_z_beta01,
        e_prior_dist01,
        e_prior_scale01,
        e_prior_df01,
        e_global_prior_df01,
        e_local01,
        e_global01,
        e_mix01,
        e_ool01,
        e_slab_df01,
        e_caux01);

beta_lp(a_z_beta01,
        a_prior_dist01,
        a_prior_scale01,
        a_prior_df01,
        a_global_prior_df01,
        a_local01,
        a_global01,
        a_mix01,
        a_ool01,
        a_slab_df01,
        a_caux01);

basehaz_lp(e_aux_unscaled01,
           e_prior_dist_for_aux01,
           e_prior_df_for_aux01);

if (e_has_intercept01 == 1)
  gamma_lp(e_gamma01[1],
           e_prior_dist_for_intercept01,
           e_prior_mean_for_intercept01,
           e_prior_scale_for_intercept01,
           e_prior_df_for_intercept01);


beta_lp(e_z_beta02,
        e_prior_dist02,
        e_prior_scale02,
        e_prior_df02,
        e_global_prior_df02,
        e_local02,
        e_global02,
        e_mix02,
        e_ool02,
        e_slab_df02,
        e_caux02);

beta_lp(a_z_beta02,
        a_prior_dist02,
        a_prior_scale02,
        a_prior_df02,
        a_global_prior_df02,
        a_local02,
        a_global02,
        a_mix02,
        a_ool02,
        a_slab_df02,
        a_caux02);

basehaz_lp(e_aux_unscaled02,
           e_prior_dist_for_aux02,
           e_prior_df_for_aux02);

if (e_has_intercept02 == 1)
  gamma_lp(e_gamma02[1],
           e_prior_dist_for_intercept02,
           e_prior_mean_for_intercept02,
           e_prior_scale_for_intercept02,
           e_prior_df_for_intercept02);


beta_lp(e_z_beta12,
        e_prior_dist12,
        e_prior_scale12,
        e_prior_df12,
        e_global_prior_df12,
        e_local12,
        e_global12,
        e_mix12,
        e_ool12,
        e_slab_df12,
        e_caux12);

beta_lp(a_z_beta12,
        a_prior_dist12,
        a_prior_scale12,
        a_prior_df12,
        a_global_prior_df12,
        a_local12,
        a_global12,
        a_mix12,
        a_ool12,
        a_slab_df12,
        a_caux12);

basehaz_lp(e_aux_unscaled12,
           e_prior_dist_for_aux12,
           e_prior_df_for_aux12);

if (e_has_intercept12 == 1)
  gamma_lp(e_gamma12[1],
           e_prior_dist_for_intercept12,
           e_prior_mean_for_intercept12,
           e_prior_scale_for_intercept12,
           e_prior_df_for_intercept12);

