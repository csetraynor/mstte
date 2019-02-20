// prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = hs
  //   4 = hs_plus
  //   5 = laplace
  //   6 = lasso
  int<lower=0,upper=6> a_prior_dist12;  // prior for assoc params

  //--- dimensions for association structure

  // num. of association parameters
  int<lower=0> a_K12;

  // used for centering assoc terms
  vector[a_K12] a_xbar12;

  // 0 = no assoc structure, 1 = any assoc structure
  int<lower=0,upper=1> assoc12;

  // which components are required to build association terms
  int<lower=0,upper=1> assoc_uses12[6,3];

  // which association terms does each submodel use
  int<lower=0,upper=1> has_assoc12[16,M];

  // num. of shared random effects
  int<lower=0> sum_size_which_b12;

  // num. of shared random effects for each long submodel
  int<lower=0> size_which_b12[M];

  // which random effects are shared for each long submodel
  int<lower=1> which_b_zindex[sum_size_which_b12];

  // num. of shared random effects incl fixed component
  int<lower=0> sum_size_which_coef12;

  // num. of shared random effects incl fixed component for each long submodel
  int<lower=0> size_which_coef12[M];

  // which random effects are shared incl fixed component
  int<lower=1> which_coef_zindex12[sum_size_which_coef12];

  // which fixed effects are shared
  int<lower=1> which_coef_xindex12[sum_size_which_coef12];

  // total num pars used in assoc*assoc interactions
  int<lower=0> sum_size_which_interactions12;

  // num pars used in assoc*assoc interactions, by submodel
  //   and by evev/evmv/mvev/mvmv interactions
  int<lower=0,upper=sum_size_which_interactions12> size_which_interactions12[M*4];

  // which terms to interact with
  int<lower=1> which_interactions12[sum_size_which_interactions12];

  // num. rows in long. predictor matrix
  int<lower=0> y_qrows12[3];

  //---- data for calculating eta in GK quadrature

  // fe design matrices

  matrix[assoc_uses12[1,1] == 1 ? y_qrows12[1] : 0, yK12[1]] y1_x_eta12;
  matrix[assoc_uses12[1,2] == 1 ? y_qrows12[2] : 0, yK12[2]] y2_x_eta12;
  matrix[assoc_uses12[1,3] == 1 ? y_qrows12[3] : 0, yK12[3]] y3_x_eta12;

  // re design matrices, group factor 1

  vector[assoc_uses12[1,1] == 1 && bK1_len[1] > 0 ? y_qrows12[1] : 0] y1_z1_eta12[bK1_len[1]];
  vector[assoc_uses12[1,2] == 1 && bK1_len[2] > 0 ? y_qrows12[2] : 0] y2_z1_eta12[bK1_len[2]];
  vector[assoc_uses12[1,3] == 1 && bK1_len[3] > 0 ? y_qrows12[3] : 0] y3_z1_eta12[bK1_len[3]];

  // re design matrices, group factor 2

  vector[assoc_uses12[1,1] == 1 && bK2_len[1] > 0 ? y_qrows12[1] : 0] y1_z2_eta12[bK2_len[1]];
  vector[assoc_uses12[1,2] == 1 && bK2_len[2] > 0 ? y_qrows12[2] : 0] y2_z2_eta12[bK2_len[2]];
  vector[assoc_uses12[1,3] == 1 && bK2_len[3] > 0 ? y_qrows12[3] : 0] y3_z2_eta12[bK2_len[3]];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z1_id_eta12[assoc_uses12[1,1] == 1 && bK1_len[1] > 0 ? y_qrows12[1] : 0];
  int<lower=0> y2_z1_id_eta12[assoc_uses12[1,2] == 1 && bK1_len[2] > 0 ? y_qrows12[2] : 0];
  int<lower=0> y3_z1_id_eta12[assoc_uses12[1,3] == 1 && bK1_len[3] > 0 ? y_qrows12[3] : 0];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z2_id_eta12[assoc_uses12[1,1] == 1 && bK2_len[1] > 0 ? y_qrows12[1] : 0];
  int<lower=0> y2_z2_id_eta12[assoc_uses12[1,2] == 1 && bK2_len[2] > 0 ? y_qrows12[2] : 0];
  int<lower=0> y3_z2_id_eta12[assoc_uses12[1,3] == 1 && bK2_len[3] > 0 ? y_qrows12[3] : 0];

  //---- data for calculating derivative of eta in GK quadrature

  // fe design matrices

  matrix[assoc_uses12[2,1] == 1 ? y_qrows12[1] : 0, yK[1]] y1_x_eps_12;
  matrix[assoc_uses12[2,2] == 1 ? y_qrows12[2] : 0, yK[2]] y2_x_eps_12;
  matrix[assoc_uses12[2,3] == 1 ? y_qrows12[3] : 0, yK[3]] y3_x_eps_12;

  // re design matrices, group factor 1

  vector[assoc_uses12[2,1] == 1 && bK1_len[1] > 0 ? y_qrows12[1] : 0] y1_z1_eps_12[bK1_len[1]];
  vector[assoc_uses12[2,2] == 1 && bK1_len[2] > 0 ? y_qrows12[2] : 0] y2_z1_eps_12[bK1_len[2]];
  vector[assoc_uses12[2,3] == 1 && bK1_len[3] > 0 ? y_qrows12[3] : 0] y3_z1_eps_12[bK1_len[3]];

  // re design matrices, group factor 2

  vector[assoc_uses12[2,1] == 1 && bK2_len[1] > 0 ? y_qrows12[1] : 0] y1_z2_eps_12[bK2_len[1]];
  vector[assoc_uses12[2,2] == 1 && bK2_len[2] > 0 ? y_qrows12[2] : 0] y2_z2_eps_12[bK2_len[2]];
  vector[assoc_uses12[2,3] == 1 && bK2_len[3] > 0 ? y_qrows12[3] : 0] y3_z2_eps_12[bK2_len[3]];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z1_id_eps_12[assoc_uses12[2,1] == 1 && bK1_len[1] > 0 ? y_qrows12[1] : 0];
  int<lower=0> y2_z1_id_eps_12[assoc_uses12[2,2] == 1 && bK1_len[2] > 0 ? y_qrows12[2] : 0];
  int<lower=0> y3_z1_id_eps_12[assoc_uses12[2,3] == 1 && bK1_len[3] > 0 ? y_qrows12[3] : 0];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z2_id_eps_12[assoc_uses12[2,1] == 1 && bK2_len[1] > 0 ? y_qrows12[1] : 0];
  int<lower=0> y2_z2_id_eps_12[assoc_uses12[2,2] == 1 && bK2_len[2] > 0 ? y_qrows12[2] : 0];
  int<lower=0> y3_z2_id_eps_12[assoc_uses12[2,3] == 1 && bK2_len[3] > 0 ? y_qrows12[3] : 0];

  //---- data for calculating integral of eta in GK quadrature

  int<lower=0> auc_qnodes12;      // num. of quadnodes for AUC of biomarker trajectory
  int<lower=0> y_qrows_for_auc12; // num. rows in long. predictor matrix at auc qpts
  vector[sum(assoc_uses12[3,]) > 0 ? y_qrows_for_auc12 : 0] auc_qwts;

  // fe design matrices

  matrix[assoc_uses12[3,1] == 1 ? y_qrows_for_auc12 : 0, yK[1]] y1_x_auc12;
  matrix[assoc_uses12[3,2] == 1 ? y_qrows_for_auc12 : 0, yK[2]] y2_x_auc12;
  matrix[assoc_uses12[3,3] == 1 ? y_qrows_for_auc12 : 0, yK[3]] y3_x_auc12;

  // re design matrices, group factor 1

  vector[assoc_uses12[3,1] == 1 && bK1_len[1] > 0 ? y_qrows_for_auc12 : 0] y1_z1_auc_12[bK1_len[1]];
  vector[assoc_uses12[3,2] == 1 && bK1_len[2] > 0 ? y_qrows_for_auc12 : 0] y2_z1_auc_12[bK1_len[2]];
  vector[assoc_uses12[3,3] == 1 && bK1_len[3] > 0 ? y_qrows_for_auc12 : 0] y3_z1_auc_12[bK1_len[3]];

  // re design matrices, group factor 2

  vector[assoc_uses12[3,1] == 1 && bK2_len[1] > 0 ? y_qrows_for_auc12 : 0] y1_z2_auc_12[bK2_len[1]];
  vector[assoc_uses12[3,2] == 1 && bK2_len[2] > 0 ? y_qrows_for_auc12 : 0] y2_z2_auc_12[bK2_len[2]];
  vector[assoc_uses12[3,3] == 1 && bK2_len[3] > 0 ? y_qrows_for_auc12 : 0] y3_z2_auc_12[bK2_len[3]];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z1_id_auc_12[assoc_uses12[3,1] == 1 && bK1_len[1] > 0 ? y_qrows_for_auc12 : 0];
  int<lower=0> y2_z1_id_auc_12[assoc_uses12[3,2] == 1 && bK1_len[2] > 0 ? y_qrows_for_auc12 : 0];
  int<lower=0> y3_z1_id_auc_12[assoc_uses12[3,3] == 1 && bK1_len[3] > 0 ? y_qrows_for_auc12 : 0];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z2_id_auc_12[assoc_uses12[3,1] == 1 && bK2_len[1] > 0 ? y_qrows_for_auc12 : 0];
  int<lower=0> y2_z2_id_auc_12[assoc_uses12[3,2] == 1 && bK2_len[2] > 0 ? y_qrows_for_auc12 : 0];
  int<lower=0> y3_z2_id_auc_-1[assoc_uses12[3,3] == 1 && bK2_len[3] > 0 ? y_qrows_for_auc12 : 0];

  //---- data for calculating assoc*data interactions in GK quadrature

  // num assoc pars used in {ev/es/mv/ms}*data interactions
  int<lower=0,upper=a_K> a_K_data12[M*4];

  // design matrix for interacting with ev/es/mv/ms at quadpoints
  matrix[sum(y_qrows12[1:M]), sum(a_K_data12)] y_x_data12;

  // indexing specifying rows of y_x_data that correspond to each submodel
  int<lower=0> idx_data12[3,2];

  //---- data for combining lower level units clustered within patients

  int<lower=0,upper=1> has_grp12[M]; // 1 = has clustering below patient level
  int<lower=0,upper=4> grp_assoc12;  // 1 = sum, 2 = mean, 3 = min, 4 = max
  int<lower=0> idx_grp12[len_cpts12,2];

