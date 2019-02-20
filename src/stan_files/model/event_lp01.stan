vector[len_epts01] e_eta_epts01;
vector[len_qpts01] e_eta_qpts01;
vector[len_ipts01] e_eta_ipts01;

real lhaz_epts01 = 0;  // summation of log hazard at event times
real lsur_qpts01 = 0;  // summation of log surv based on qpts01
real lsur_ipts01 = 0;  // summation of log surv based on ipts01

if (len_epts01 > 0) e_eta_epts01 = e_eta01[idx_cpts01[1,1]:idx_cpts01[1,2]];
if (len_qpts01 > 0) e_eta_qpts01 = e_eta01[idx_cpts01[2,1]:idx_cpts01[2,2]];
if (len_ipts01 > 0) e_eta_ipts01 = e_eta01[idx_cpts01[3,1]:idx_cpts01[3,2]];

// evaluate log hazard and log survival
if (basehaz_type01 == 5) { // exponential model
  if (len_epts01 > 0) {
    lhaz_epts01 = sum(e_eta_epts01);
  }
  if (len_qpts01 > 0) {
    lsur_qpts01 = - dot_product(qwts01, exp(e_eta_qpts01));
  }
  if (len_ipts01 > 0) {
    lsur_ipts01 = - dot_product(iwts01, exp(e_eta_ipts01));
  }
}
else if (basehaz_type01 == 1) { // weibull model
  real shape = e_aux01[1];
  real log_shape = log(shape);
  if (len_epts01 > 0) {
    lhaz_epts01 = (len_epts01 * log_shape) + (shape - 1) * sum_log_epts01 + sum(e_eta_epts01);
  }
  if (len_qpts01 > 0) {
    vector[len_qpts01] lhaz_qpts01;
    lhaz_qpts01 = log_shape + (shape - 1) * log_qpts01 + e_eta_qpts01;
    lsur_qpts01 = - dot_product(qwts01, exp(lhaz_qpts01));
  }
  if (len_ipts01 > 0) {
    vector[len_ipts01] lhaz_ipts01;
    lhaz_ipts01 = log_shape + (shape - 1) * log_ipts01 + e_eta_ipts01;
    lsur_ipts01 = - dot_product(iwts01, exp(lhaz_ipts01));
  }
}
else if (basehaz_type01 == 6) { // gompertz model
  real scale = e_aux01[1];
  if (len_epts01 > 0) {
    lhaz_epts01 = scale * sum_epts01 + sum(e_eta_epts01);
  }
  if (len_qpts01 > 0) {
    vector[len_qpts01] lhaz_qpts01;
    lhaz_qpts01 = scale * qpts01 + e_eta_qpts01;
    lsur_qpts01 = - dot_product(qwts01, exp(lhaz_qpts01));
  }
  if (len_ipts01 > 0) {
    vector[len_ipts01] lhaz_ipts01;
    lhaz_ipts01 = scale * ipts01 + e_eta_ipts01;
    lsur_ipts01 = - dot_product(iwts01, exp(lhaz_ipts01));
  }
}
else if (basehaz_type01 == 4) { // M-splines, on haz scale
  if (len_epts01 > 0)
    lhaz_epts01 = sum(log(basis_epts01 * e_aux01) + e_eta_epts01);
  if (len_qpts01 > 0) {
    vector[len_qpts01] lhaz_qpts01;
    lhaz_qpts01 = log(basis_qpts01 * e_aux01) + e_eta_qpts01;
    lsur_qpts01 = - dot_product(qwts01, exp(lhaz_qpts01));
  }
  if (len_ipts01 > 0) {
    vector[len_ipts01] lhaz_ipts01;
    lhaz_ipts01 = log(basis_ipts01 * e_aux01) + e_eta_ipts01;
    lsur_ipts01 = - dot_product(iwts01, exp(lhaz_ipts01));
  }
}
else if (basehaz_type01 == 2) { // B-splines, on log haz scale
  if (len_epts01 > 0) {
    lhaz_epts01 = sum(basis_epts01 * e_aux01 + e_eta_epts01);
  }
  if (len_qpts01 > 0) {
    vector[len_qpts01] lhaz_qpts01;
    lhaz_qpts01 = basis_qpts01 * e_aux01 + e_eta_qpts01;
    lsur_qpts01 = - dot_product(qwts01, exp(lhaz_qpts01));
  }
  if (len_ipts01 > 0) {
    vector[len_ipts01] lhaz_ipts01;
    lhaz_ipts01 = basis_ipts01 * e_aux01 + e_eta_ipts01;
    lsur_ipts01 = - dot_product(iwts01, exp(lhaz_ipts01));
  }
}
else {
  reject("Bug found: invalid baseline hazard.");
}

// increment target with log-lik for event submodel
if (has_weights01 == 0 && prior_PD == 0) { // unweighted log likelihood
  target += lhaz_epts01 + lsur_qpts01 - lsur_ipts01;
}
else if (prior_PD == 0) { // weighted log likelihood
  reject("Bug found: weights are not yet implemented.");
}
