vector[len_epts12] e_eta_epts12;
vector[len_qpts12] e_eta_qpts12;
vector[len_ipts12] e_eta_ipts12;

real lhaz_epts12 = 0;  // summation of log hazard at event times
real lsur_qpts12 = 0;  // summation of log surv based on qpts12
real lsur_ipts12 = 0;  // summation of log surv based on ipts12

if (len_epts12 > 0) e_eta_epts12 = e_eta12[idx_cpts12[1,1]:idx_cpts12[1,2]];
if (len_qpts12 > 0) e_eta_qpts12 = e_eta12[idx_cpts12[2,1]:idx_cpts12[2,2]];
if (len_ipts12 > 0) e_eta_ipts12 = e_eta12[idx_cpts12[3,1]:idx_cpts12[3,2]];

// evaluate log hazard and log survival
if (basehaz_type12 == 5) { // exponential model
  if (len_epts12 > 0) {
    lhaz_epts12 = sum(e_eta_epts12);
  }
  if (len_qpts12 > 0) {
    lsur_qpts12 = - dot_product(qwts12, exp(e_eta_qpts12));
  }
  if (len_ipts12 > 0) {
    lsur_ipts12 = - dot_product(iwts12, exp(e_eta_ipts12));
  }
}
else if (basehaz_type12 == 1) { // weibull model
  real shape = e_aux12[1];
  real log_shape = log(shape);
  if (len_epts12 > 0) {
    lhaz_epts12 = (len_epts12 * log_shape) + (shape - 1) * sum_log_epts12 + sum(e_eta_epts12);
  }
  if (len_qpts12 > 0) {
    vector[len_qpts12] lhaz_qpts12;
    lhaz_qpts12 = log_shape + (shape - 1) * log_qpts12 + e_eta_qpts12;
    lsur_qpts12 = - dot_product(qwts12, exp(lhaz_qpts12));
  }
  if (len_ipts12 > 0) {
    vector[len_ipts12] lhaz_ipts12;
    lhaz_ipts12 = log_shape + (shape - 1) * log_ipts12 + e_eta_ipts12;
    lsur_ipts12 = - dot_product(iwts12, exp(lhaz_ipts12));
  }
}
else if (basehaz_type12 == 6) { // gompertz model
  real scale = e_aux12[1];
  if (len_epts12 > 0) {
    lhaz_epts12 = scale * sum_epts12 + sum(e_eta_epts12);
  }
  if (len_qpts12 > 0) {
    vector[len_qpts12] lhaz_qpts12;
    lhaz_qpts12 = scale * qpts12 + e_eta_qpts12;
    lsur_qpts12 = - dot_product(qwts12, exp(lhaz_qpts12));
  }
  if (len_ipts12 > 0) {
    vector[len_ipts12] lhaz_ipts12;
    lhaz_ipts12 = scale * ipts12 + e_eta_ipts12;
    lsur_ipts12 = - dot_product(iwts12, exp(lhaz_ipts12));
  }
}
else if (basehaz_type12 == 4) { // M-splines, on haz scale
  if (len_epts12 > 0)
    lhaz_epts12 = sum(log(basis_epts12 * e_aux12) + e_eta_epts12);
  if (len_qpts12 > 0) {
    vector[len_qpts12] lhaz_qpts12;
    lhaz_qpts12 = log(basis_qpts12 * e_aux12) + e_eta_qpts12;
    lsur_qpts12 = - dot_product(qwts12, exp(lhaz_qpts12));
  }
  if (len_ipts12 > 0) {
    vector[len_ipts12] lhaz_ipts12;
    lhaz_ipts12 = log(basis_ipts12 * e_aux12) + e_eta_ipts12;
    lsur_ipts12 = - dot_product(iwts12, exp(lhaz_ipts12));
  }
}
else if (basehaz_type12 == 2) { // B-splines, on log haz scale
  if (len_epts12 > 0) {
    lhaz_epts12 = sum(basis_epts12 * e_aux12 + e_eta_epts12);
  }
  if (len_qpts12 > 0) {
    vector[len_qpts12] lhaz_qpts12;
    lhaz_qpts12 = basis_qpts12 * e_aux12 + e_eta_qpts12;
    lsur_qpts12 = - dot_product(qwts12, exp(lhaz_qpts12));
  }
  if (len_ipts12 > 0) {
    vector[len_ipts12] lhaz_ipts12;
    lhaz_ipts12 = basis_ipts12 * e_aux12 + e_eta_ipts12;
    lsur_ipts12 = - dot_product(iwts12, exp(lhaz_ipts12));
  }
}
else {
  reject("Bug found: invalid baseline hazard.");
}

// increment target with log-lik for event submodel
if (has_weights == 0 && prior_PD == 0) { // unweighted log likelihood
  target += lhaz_epts12 + lsur_qpts12 - lsur_ipts12;
}
else if (prior_PD == 0) { // weighted log likelihood
  reject("Bug found: weights are not yet implemented.");
}
