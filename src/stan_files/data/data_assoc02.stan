// prior family:
  //   0 = none
  //   1 = normal
  //   2 = student_t
  //   3 = hs
  //   4 = hs_plus
  //   5 = laplace
  //   6 = lasso
  int<lower=0,upper=6> a_prior_dist02;  // prior for assoc params

  //--- dimensions for association structure

  // num. of association parameters
  int<lower=0> a_K02;

  // used for centering assoc terms
  vector[a_K02] a_xbar02;

  // 0 = no assoc structure, 1 = any assoc structure
  int<lower=0,upper=1> assoc02;

  // which components are required to build association terms
  int<lower=0,upper=1> assoc_uses02[6,3];

  // which association terms does each submodel use
  int<lower=0,upper=1> has_assoc02[16,M];

  // num. of shared random effects
  int<lower=0> sum_size_which_b02;

  // num. of shared random effects for each long submodel
  int<lower=0> size_which_b02[M];

  // which random effects are shared for each long submodel
  int<lower=1> which_b_zindex02[sum_size_which_b02];

  // num. of shared random effects incl fixed component
  int<lower=0> sum_size_which_coef02;

  // num. of shared random effects incl fixed component for each long submodel
  int<lower=0> size_which_coef02[M];

  // which random effects are shared incl fixed component
  int<lower=1> which_coef_zindex02[sum_size_which_coef02];

  // which fixed effects are shared
  int<lower=1> which_coef_xindex02[sum_size_which_coef02];

  // total num pars used in assoc*assoc interactions
  int<lower=0> sum_size_which_interactions02;

  // num pars used in assoc*assoc interactions, by submodel
  //   and by evev/evmv/mvev/mvmv interactions
  int<lower=0,upper=sum_size_which_interactions02> size_which_interactions02[M*4];

  // which terms to interact with
  int<lower=1> which_interactions02[sum_size_which_interactions02];

  // num. rows in long. predictor matrix
  int<lower=0> y_qrows02[3];

  //---- data for calculating eta in GK quadrature

  // fe design matrices

  matrix[assoc_uses02[1,1] == 1 ? y_qrows02[1] : 0, yK[1]] y1_x_eta_02;
  matrix[assoc_uses02[1,2] == 1 ? y_qrows02[2] : 0, yK[2]] y2_x_eta_02;
  matrix[assoc_uses02[1,3] == 1 ? y_qrows02[3] : 0, yK[3]] y3_x_eta_02;

  // re design matrices, group factor 1

  vector[assoc_uses02[1,1] == 1 && bK1_len[1] > 0 ? y_qrows02[1] : 0] y1_z1_eta_02[bK1_len[1]];
  vector[assoc_uses02[1,2] == 1 && bK1_len[2] > 0 ? y_qrows02[2] : 0] y2_z1_eta_02[bK1_len[2]];
  vector[assoc_uses02[1,3] == 1 && bK1_len[3] > 0 ? y_qrows02[3] : 0] y3_z1_eta_02[bK1_len[3]];

  // re design matrices, group factor 2

  vector[assoc_uses02[1,1] == 1 && bK2_len[1] > 0 ? y_qrows02[1] : 0] y1_z2_eta_02[bK2_len[1]];
  vector[assoc_uses02[1,2] == 1 && bK2_len[2] > 0 ? y_qrows02[2] : 0] y2_z2_eta_02[bK2_len[2]];
  vector[assoc_uses02[1,3] == 1 && bK2_len[3] > 0 ? y_qrows02[3] : 0] y3_z2_eta_02[bK2_len[3]];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z1_id_eta_02[assoc_uses02[1,1] == 1 && bK1_len[1] > 0 ? y_qrows02[1] : 0];
  int<lower=0> y2_z1_id_eta_02[assoc_uses02[1,2] == 1 && bK1_len[2] > 0 ? y_qrows02[2] : 0];
  int<lower=0> y3_z1_id_eta_02[assoc_uses02[1,3] == 1 && bK1_len[3] > 0 ? y_qrows02[3] : 0];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z2_id_eta_02[assoc_uses02[1,1] == 1 && bK2_len[1] > 0 ? y_qrows02[1] : 0];
  int<lower=0> y2_z2_id_eta_02[assoc_uses02[1,2] == 1 && bK2_len[2] > 0 ? y_qrows02[2] : 0];
  int<lower=0> y3_z2_id_eta_02[assoc_uses02[1,3] == 1 && bK2_len[3] > 0 ? y_qrows02[3] : 0];

  //---- data for calculating derivative of eta in GK quadrature

  // fe design matrices

  matrix[assoc_uses02[2,1] == 1 ? y_qrows02[1] : 0, yK[1]] y1_x_eps_02;
  matrix[assoc_uses02[2,2] == 1 ? y_qrows02[2] : 0, yK[2]] y2_x_eps_02;
  matrix[assoc_uses02[2,3] == 1 ? y_qrows02[3] : 0, yK[3]] y3_x_eps_02;

  // re design matrices, group factor 1

  vector[assoc_uses02[2,1] == 1 && bK1_len[1] > 0 ? y_qrows02[1] : 0] y1_z1_eps_02[bK1_len[1]];
  vector[assoc_uses02[2,2] == 1 && bK1_len[2] > 0 ? y_qrows02[2] : 0] y2_z1_eps_02[bK1_len[2]];
  vector[assoc_uses02[2,3] == 1 && bK1_len[3] > 0 ? y_qrows02[3] : 0] y3_z1_eps_02[bK1_len[3]];

  // re design matrices, group factor 2

  vector[assoc_uses02[2,1] == 1 && bK2_len[1] > 0 ? y_qrows02[1] : 0] y1_z2_eps_02[bK2_len[1]];
  vector[assoc_uses02[2,2] == 1 && bK2_len[2] > 0 ? y_qrows02[2] : 0] y2_z2_eps_02[bK2_len[2]];
  vector[assoc_uses02[2,3] == 1 && bK2_len[3] > 0 ? y_qrows02[3] : 0] y3_z2_eps_02[bK2_len[3]];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z1_id_eps_02[assoc_uses02[2,1] == 1 && bK1_len[1] > 0 ? y_qrows02[1] : 0];
  int<lower=0> y2_z1_id_eps_02[assoc_uses02[2,2] == 1 && bK1_len[2] > 0 ? y_qrows02[2] : 0];
  int<lower=0> y3_z1_id_eps_02[assoc_uses02[2,3] == 1 && bK1_len[3] > 0 ? y_qrows02[3] : 0];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z2_id_eps_02[assoc_uses02[2,1] == 1 && bK2_len[1] > 0 ? y_qrows02[1] : 0];
  int<lower=0> y2_z2_id_eps_02[assoc_uses02[2,2] == 1 && bK2_len[2] > 0 ? y_qrows02[2] : 0];
  int<lower=0> y3_z2_id_eps_02[assoc_uses02[2,3] == 1 && bK2_len[3] > 0 ? y_qrows02[3] : 0];

  //---- data for calculating integral of eta in GK quadrature

  int<lower=0> auc_qnodes02;      // num. of quadnodes for AUC of biomarker trajectory
  int<lower=0> y_qrows_for_auc02; // num. rows in long. predictor matrix at auc qpts
  vector[sum(assoc_uses02[3,]) > 0 ? y_qrows_for_auc02 : 0] auc_qwts02;

  // fe design matrices

  matrix[assoc_uses02[3,1] == 1 ? y_qrows_for_auc02 : 0, yK[1]] y1_x_auc_02;
  matrix[assoc_uses02[3,2] == 1 ? y_qrows_for_auc02 : 0, yK[2]] y2_x_auc_02;
  matrix[assoc_uses02[3,3] == 1 ? y_qrows_for_auc02 : 0, yK[3]] y3_x_auc_02;

  // re design matrices, group factor 1

  vector[assoc_uses02[3,1] == 1 && bK1_len[1] > 0 ? y_qrows_for_auc02 : 0] y1_z1_auc_02[bK1_len[1]];
  vector[assoc_uses02[3,2] == 1 && bK1_len[2] > 0 ? y_qrows_for_auc02 : 0] y2_z1_auc_02[bK1_len[2]];
  vector[assoc_uses02[3,3] == 1 && bK1_len[3] > 0 ? y_qrows_for_auc02 : 0] y3_z1_auc_02[bK1_len[3]];

  // re design matrices, group factor 2

  vector[assoc_uses02[3,1] == 1 && bK2_len[1] > 0 ? y_qrows_for_auc02 : 0] y1_z2_auc_02[bK2_len[1]];
  vector[assoc_uses02[3,2] == 1 && bK2_len[2] > 0 ? y_qrows_for_auc02 : 0] y2_z2_auc_02[bK2_len[2]];
  vector[assoc_uses02[3,3] == 1 && bK2_len[3] > 0 ? y_qrows_for_auc02 : 0] y3_z2_auc_02[bK2_len[3]];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z1_id_auc_02[assoc_uses02[3,1] == 1 && bK1_len[1] > 0 ? y_qrows_for_auc02 : 0];
  int<lower=0> y2_z1_id_auc_02[assoc_uses02[3,2] == 1 && bK1_len[2] > 0 ? y_qrows_for_auc02 : 0];
  int<lower=0> y3_z1_id_auc_02[assoc_uses02[3,3] == 1 && bK1_len[3] > 0 ? y_qrows_for_auc02 : 0];

  // ids for re design matrices, group factor 1

  int<lower=0> y1_z2_id_auc_02[assoc_uses02[3,1] == 1 && bK2_len[1] > 0 ? y_qrows_for_auc02 : 0];
  int<lower=0> y2_z2_id_auc_02[assoc_uses02[3,2] == 1 && bK2_len[2] > 0 ? y_qrows_for_auc02 : 0];
  int<lower=0> y3_z2_id_auc_02[assoc_uses02[3,3] == 1 && bK2_len[3] > 0 ? y_qrows_for_auc02 : 0];

  //---- data for calculating assoc*data interactions in GK quadrature

  // num assoc pars used in {ev/es/mv/ms}*data interactions
  int<lower=0,upper=a_K02> a_K_data02[M*4];

  // design matrix for interacting with ev/es/mv/ms at quadpoints
  matrix[sum(y_qrows02[1:M]), sum(a_K_data02)] y_x_data02;

  // indexing specifying rows of y_x_data that correspond to each submodel
  int<lower=0> idx_data02[3,2];

  //---- data for combining lower level units clustered within patients

  int<lower=0,upper=1> has_grp02[M]; // 1 = has clustering below patient level
  int<lower=0,upper=4> grp_assoc02;  // 1 = sum, 2 = mean, 3 = min, 4 = max
  int<lower=0> idx_grp02[len_cpts02,2];

