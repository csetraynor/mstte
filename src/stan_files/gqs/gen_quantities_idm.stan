// norm_const is a constant shift in log baseline hazard
if (e_has_intercept01 == 1)
  e_alpha01 = e_gamma01[1] + norm_const01 -
    dot_product(e_xbar01, e_beta01) - dot_product(a_xbar01, a_beta01);
else
  e_alpha01 = norm_const01 -
    dot_product(e_xbar01, e_beta01) - dot_product(a_xbar01, a_beta01);

// norm_const is a constant shift in log baseline hazard
if (e_has_intercept02 == 1)
  e_alpha02 = e_gamma02[1] + norm_const02 -
    dot_product(e_xbar02, e_beta02) - dot_product(a_xbar02, a_beta02);
else
  e_alpha02 = norm_const02 -
    dot_product(e_xbar02, e_beta02) - dot_product(a_xbar02, a_beta02);

// norm_const is a constant shift in log baseline hazard
if (e_has_intercept12 == 1)
  e_alpha12 = e_gamma12[1] + norm_const12 -
    dot_product(e_xbar12, e_beta12) - dot_product(a_xbar12, a_beta12);
else
  e_alpha12 = norm_const12 -
    dot_product(e_xbar12, e_beta12) - dot_product(a_xbar12, a_beta12);
