vector[len_epts02] e_eta_epts02;
vector[len_qpts02] e_eta_qpts02;
vector[len_ipts02] e_eta_ipts02;

real lhaz_epts02 = 0;  // summation of log hazard at event times
real lsur_qpts02 = 0;  // summation of log surv based on qpts02
real lsur_ipts02 = 0;  // summation of log surv based on ipts02

if (len_epts02 > 0) e_eta_epts02 = e_eta02[idx_cpts02[1,1]:idx_cpts02[1,2]];
if (len_qpts02 > 0) e_eta_qpts02 = e_eta02[idx_cpts02[2,1]:idx_cpts02[2,2]];
if (len_ipts02 > 0) e_eta_ipts02 = e_eta02[idx_cpts02[3,1]:idx_cpts02[3,2]];

// evaluate log hazard and log survival
if (basehaz_type02 == 5) { // exponential model
  if (len_epts02 > 0) {
    lhaz_epts02 = sum(e_eta_epts02);
  }
  if (len_qpts02 > 0) {
    lsur_qpts02 = - dot_product(qwts02, exp(e_eta_qpts02));
  }
  if (len_ipts02 > 0) {
    lsur_ipts02 = - dot_product(iwts02, exp(e_eta_ipts02));
  }
}
else if (basehaz_type02 == 1) { // weibull model
  real shape = e_aux02[1];
  real log_shape = log(shape);
  if (len_epts02 > 0) {
    lhaz_epts02 = (len_epts02 * log_shape) + (shape - 1) * sum_log_epts02 + sum(e_eta_epts02);
  }
  if (len_qpts02 > 0) {
    vector[len_qpts02] lhaz_qpts02;
    lhaz_qpts02 = log_shape + (shape - 1) * log_qpts02 + e_eta_qpts02;
    lsur_qpts02 = - dot_product(qwts02, exp(lhaz_qpts02));
  }
  if (len_ipts02 > 0) {
    vector[len_ipts02] lhaz_ipts02;
    lhaz_ipts02 = log_shape + (shape - 1) * log_ipts02 + e_eta_ipts02;
    lsur_ipts02 = - dot_product(iwts02, exp(lhaz_ipts02));
  }
}
else if (basehaz_type02 == 6) { // gompertz model
  real scale = e_aux02[1];
  if (len_epts02 > 0) {
    lhaz_epts02 = scale * sum_epts02 + sum(e_eta_epts02);
  }
  if (len_qpts02 > 0) {
    vector[len_qpts02] lhaz_qpts02;
    lhaz_qpts02 = scale * qpts02 + e_eta_qpts02;
    lsur_qpts02 = - dot_product(qwts02, exp(lhaz_qpts02));
  }
  if (len_ipts02 > 0) {
    vector[len_ipts02] lhaz_ipts02;
    lhaz_ipts02 = scale * ipts02 + e_eta_ipts02;
    lsur_ipts02 = - dot_product(iwts02, exp(lhaz_ipts02));
  }
}
else if (basehaz_type02 == 4) { // M-splines, on haz scale
  if (len_epts02 > 0)
    lhaz_epts02 = sum(log(basis_epts02 * e_aux02) + e_eta_epts02);
  if (len_qpts02 > 0) {
    vector[len_qpts02] lhaz_qpts02;
    lhaz_qpts02 = log(basis_qpts02 * e_aux02) + e_eta_qpts02;
    lsur_qpts02 = - dot_product(qwts02, exp(lhaz_qpts02));
  }
  if (len_ipts02 > 0) {
    vector[len_ipts02] lhaz_ipts02;
    lhaz_ipts02 = log(basis_ipts02 * e_aux02) + e_eta_ipts02;
    lsur_ipts02 = - dot_product(iwts02, exp(lhaz_ipts02));
  }
}
else if (basehaz_type02 == 2) { // B-splines, on log haz scale
  if (len_epts02 > 0) {
    lhaz_epts02 = sum(basis_epts02 * e_aux02 + e_eta_epts02);
  }
  if (len_qpts02 > 0) {
    vector[len_qpts02] lhaz_qpts02;
    lhaz_qpts02 = basis_qpts02 * e_aux02 + e_eta_qpts02;
    lsur_qpts02 = - dot_product(qwts02, exp(lhaz_qpts02));
  }
  if (len_ipts02 > 0) {
    vector[len_ipts02] lhaz_ipts02;
    lhaz_ipts02 = basis_ipts02 * e_aux02 + e_eta_ipts02;
    lsur_ipts02 = - dot_product(iwts02, exp(lhaz_ipts02));
  }
}
else {
  reject("Bug found: invalid baseline hazard.");
}

// increment target with log-lik for event submodel
if (has_weights == 0 && prior_PD == 0) { // unweighted log likelihood
  target += lhaz_epts02 + lsur_qpts02 - lsur_ipts02;
}
else if (prior_PD == 0) { // weighted log likelihood
  reject("Bug found: weights are not yet implemented.");
}
