    // !!! be careful that indexing of has_assoc01 matches stan_jm.fit !!!

    // mark tracks indexing within a_beta01 vector, which is the
    // vector of association parameters
    int mark = 0;

    // mark2 tracks indexing within a_K_data vector, which is the
    // vector specifying the number of columns used for each possible
    // type of association term by data interaction
    int mark2 = 0;

    // mark3 tracks indexing within size_which_interactions vector
    int mark3 = 0;

    for (m in 1:M) {

      //----- etavalue and any interactions

      mark2 += 1;
      if (has_assoc01[1,m]  == 1 || // etavalue
          has_assoc01[9,m]  == 1 || // etavalue * data
          has_assoc01[13,m] == 1 || // etavalue * etavalue
          has_assoc01[14,m] == 1) { // etavalue * muvalue

        // declare and define eta for longitudinal submodel m
#include /model/make_eta_tmp01.stan

        // add etavalue and any interactions to event submodel eta
        if (has_assoc01[1,m] == 1) { // etavalue
          vector[y_qrows01[m]] val;
          if (has_grp01[m] == 0) {
            val = eta_tmp01;
          }
          else {
            val = collapse_within_groups(eta_tmp01, idx_grp01, grp_assoc01);
          }
          mark += 1;
          e_eta01 += a_beta01[mark] * (val - a_xbar01[mark]);
        }
        // !! Note that for now I only have assured that eta value will work!
        /*

        mark2 += 1; // count even if assoc type isn't used
        if (has_assoc01[9,m] == 1) { // etavalue*data
          int idx1 = idx_data01[m,1];
          int idx2 = idx_data01[m,2];
          int J = a_K_data[mark2];
          int j_shift = (mark2 == 1) ? 0 : sum(a_K_data[1:(mark2-1)]);
          for (j in 1:J) {
            vector[y_qrows01[m]] val;
            int sel = j_shift + j;
            if (has_grp01[m] == 0) {
              val = eta_tmp .* y_x_data01[idx1:idx2, sel];
            }
            else {
              val = collapse_within_groups(
                eta_tmp .* y_x_data01[idx1:idx2, sel],
                idx_grp01, grp_assoc01);
            }
            mark += 1;
            e_eta01 += a_beta01[mark] * (val - a_xbar01[mark]);
          }
        }
        mark3 += 1; // count even if assoc type isn't used
        if (has_assoc01[13,m] == 1) { // etavalue*etavalue
          for (j in 1:size_which_interactions[mark3]) {
            int j_shift = (mark3 == 1) ? 0 : sum(size_which_interactions[1:(mark3-1)]);
            int sel = which_interactions[j+j_shift];
            vector[y_qrows01[m]] val;
#include /model/make_eta_tmp2.stan
            val = eta_tmp .* eta_tmp2;
            mark += 1;
            e_eta01 += a_beta01[mark] * (val - a_xbar01[mark]);
          }
        }
        mark3 += 1; // count even if assoc type isn't used
        if (has_assoc01[14,m] == 1) { // etavalue*muvalue
          for (j in 1:size_which_interactions[mark3]) {
            int j_shift = (mark3 == 1) ? 0 : sum(size_which_interactions[1:(mark3-1)]);
            int sel = which_interactions[j+j_shift];
            vector[y_qrows01[m]] val;
            vector[y_qrows01[sel]] mu_tmp2;
#include /model/make_eta_tmp2.stan
            mu_tmp2 = evaluate_mu(eta_tmp2, family[sel], link[sel]);
            val = eta_tmp .* mu_tmp2;
            mark += 1;
            e_eta01 += a_beta01[mark] * (val - a_xbar01[mark]);
          }
        }
      }
      else {
        mark3 += 2;
      }

      //----- etaslope and any interactions

      mark2 += 1;
      if ((has_assoc01[2,m] == 1) || (has_assoc01[10,m] == 1)) {

        // declare and define etaslope at quadpoints for submodel m
        vector[y_qrows01[m]] dydt_eta_q;
        if (m == 1) {
          int bMat1_colshift = 0;
          int bMat2_colshift = 0;
          dydt_eta_q = evaluate_eta01(y1_x_eps_01,
                                    y1_z1_eps_01,
                                    y1_z2_eps_01,
                                    y1_z1_id_eps_01,
                                    y1_z2_id_eps_01,
                                    yGamma1,
                                    yBeta1,
                                    bMat1,
                                    bMat2,
                                    bMat1_colshift,
                                    bMat2_colshift,
                                    0);
        }
        else if (m == 2) {
          int bMat1_colshift = bK1_len[1];
          int bMat2_colshift = bK2_len[1];
          dydt_eta_q = evaluate_eta01(y2_x_eps_01,
                                    y2_z1_eps_01,
                                    y2_z2_eps_01,
                                    y2_z1_id_eps_01,
                                    y2_z2_id_eps_01,
                                    yGamma2,
                                    yBeta2,
                                    bMat1,
                                    bMat2,
                                    bMat1_colshift,
                                    bMat2_colshift,
                                    0);
        }
        else if (m == 3) {
          int bMat1_colshift = sum(bK1_len[1:2]);
          int bMat2_colshift = sum(bK2_len[1:2]);
          dydt_eta_q = evaluate_eta01(y3_x_eps_01,
                                    y3_z1_eps_01,
                                    y3_z2_eps_01,
                                    y3_z1_id_eps_01,
                                    y3_z2_id_eps_01,
                                    yGamma3,
                                    yBeta3,
                                    bMat1,
                                    bMat2,
                                    bMat1_colshift,
                                    bMat2_colshift,
                                    0);
        }

        // add etaslope and any interactions to event submodel eta
        if (has_assoc01[2,m] == 1) { // etaslope
          vector[y_qrows01[m]] val;
          if (has_grp01[m] == 0) {
            val = dydt_eta_q;
          }
          else {
            val = collapse_within_groups(dydt_eta_q, idx_grp01, grp_assoc01);
          }
          mark += 1;
          e_eta01 = e_eta01 + a_beta01[mark] * (val - a_xbar01[mark]);
        }
        if (has_assoc01[10,m] == 1) { // etaslope*data
          int idx1 = idx_data01[m,1];
          int idx2 = idx_data01[m,2];
          int J = a_K_data[mark2];
          int j_shift = (mark2 == 1) ? 0 : sum(a_K_data[1:(mark2-1)]);
          for (j in 1:J) {
            vector[y_qrows01[m]] val;
            int sel = j_shift + j;
            if (has_grp01[m] == 0) {
              val = dydt_eta_q .* y_x_data01[idx1:idx2, sel];
            }
            else {
              val = collapse_within_groups(dydt_eta_q .* y_x_data01[idx1:idx2, sel],
                                           idx_grp01, grp_assoc01);
            }
            mark += 1;
            e_eta01 = e_eta01 + a_beta01[mark] * (val - a_xbar01[mark]);
          }
        }
      }

      //----- etaauc

      if (has_assoc01[3,m] == 1) { // etaauc
        vector[y_qrows_for_auc01] eta_auc_tmp; // eta at all auc qpts (for submodel m)
        vector[y_qrows01[m]] val;              // eta following summation over auc qpts
        if (m == 1) {
          int bMat1_colshift = 0;
          int bMat2_colshift = 0;
          eta_auc_tmp = evaluate_eta01(y1_x_auc_01,
                                     y1_z1_auc_01,
                                     y1_z2_auc_01,
                                     y1_z1_id_auc_01,
                                     y1_z2_id_auc_01,
                                     yGamma1,
                                     yBeta1,
                                     bMat1,
                                     bMat2,
                                     bMat1_colshift,
                                     bMat2_colshift,
                                     intercept_type[1]);
        }
        else if (m == 2) {
          int bMat1_colshift = bK1_len[1];
          int bMat2_colshift = bK2_len[1];
          eta_auc_tmp = evaluate_eta01(y2_x_auc_01,
                                     y2_z1_auc_01,
                                     y2_z2_auc_01,
                                     y2_z1_id_auc_01,
                                     y2_z2_id_auc_01,
                                     yGamma2,
                                     yBeta2,
                                     bMat1,
                                     bMat2,
                                     bMat1_colshift,
                                     bMat2_colshift,
                                     intercept_type[2]);
        }
        else if (m == 3) {
          int bMat1_colshift = sum(bK1_len[1:2]);
          int bMat2_colshift = sum(bK2_len[1:2]);
          eta_auc_tmp = evaluate_eta01(y3_x_auc_01,
                                     y3_z1_auc_01,
                                     y3_z2_auc_01,
                                     y3_z1_id_auc_01,
                                     y3_z2_id_auc_01,
                                     yGamma3,
                                     yBeta3,
                                     bMat1,
                                     bMat2,
                                     bMat1_colshift,
                                     bMat2_colshift,
                                     intercept_type[3]);
        }
        mark += 1;
        for (r in 1:y_qrows01[m]) {
          vector[auc_qnodes01] val_tmp;
          vector[auc_qnodes01] wgt_tmp;
          val_tmp = eta_auc_tmp[((r-1) * auc_qnodes01 + 1):(r * auc_qnodes01)];
          wgt_tmp = auc_qwts   [((r-1) * auc_qnodes01 + 1):(r * auc_qnodes01)];
          val[r] = sum(wgt_tmp .* val_tmp);
        }
        e_eta01 += a_beta01[mark] * (val - a_xbar01[mark]);
      }

      //----- muvalue and any interactions

      mark2 += 1;
      if (has_assoc01[4,m]  == 1 || // muvalue
          has_assoc01[11,m] == 1 || // muvalue * data
          has_assoc01[15,m] == 1 || // muvalue * etavalue
          has_assoc01[16,m] == 1) { // muvalue * muvalue

        // declare and define mu for submodel m
        vector[y_qrows01[m]] mu_tmp;
#include /model/make_eta_tmp.stan
        mu_tmp = evaluate_mu(eta_tmp, family[m], link[m]);

        // add muvalue and any interactions to event submodel eta
        if (has_assoc01[4,m] == 1) { // muvalue
          vector[y_qrows01[m]] val;
          if (has_grp01[m] == 0) {
            val = mu_tmp;
          }
          else {
            val = collapse_within_groups(mu_tmp, idx_grp01, grp_assoc01);
          }
          mark += 1;
          e_eta01 = e_eta01 + a_beta01[mark] * (val - a_xbar01[mark]);
        }
        if (has_assoc01[11,m] == 1) { // muvalue*data
          int tmp = a_K_data[mark2];
          int j_shift = (mark2 == 1) ? 0 : sum(a_K_data[1:(mark2-1)]);
          for (j in 1:tmp) {
            vector[y_qrows01[m]] val;
            int sel = j_shift + j;
            if (has_grp01[m] == 0) {
              val = mu_tmp .* y_x_data01[idx_data01[m,1]:idx_data01[m,2], sel];
            }
            else {
              val = collapse_within_groups(
                mu_tmp .* y_x_data01[idx_data01[m,1]:idx_data01[m,2], sel],
                idx_grp01, grp_assoc01);
            }
            mark += 1;
            e_eta01 = e_eta01 + a_beta01[mark] * (val - a_xbar01[mark]);
          }
        }
        mark3 += 1; // count even if assoc type isn't used
        if (has_assoc01[15,m] == 1) { // muvalue*etavalue
          for (j in 1:size_which_interactions[mark3]) {
            int j_shift = (mark3 == 1) ? 0 : sum(size_which_interactions[1:(mark3-1)]);
            int sel = which_interactions[j+j_shift];
            vector[y_qrows01[m]] val;
#include /model/make_eta_tmp2.stan
            val = mu_tmp .* eta_tmp2;
            mark += 1;
            e_eta01 = e_eta01 + a_beta01[mark] * (val - a_xbar01[mark]);
          }
        }
        mark3 += 1; // count even if assoc type isn't used
        if (has_assoc01[16,m] == 1) { // muvalue*muvalue
          for (j in 1:size_which_interactions[mark3]) {
            int j_shift = (mark3 == 1) ? 0 : sum(size_which_interactions[1:(mark3-1)]);
            int sel = which_interactions[j+j_shift];
            vector[y_qrows01[m]] val;
            vector[y_qrows01[sel]] mu_tmp2;
#include /model/make_eta_tmp2.stan
            mu_tmp2 = evaluate_mu(eta_tmp2, family[sel], link[sel]);
            val = mu_tmp .* mu_tmp2;
            mark += 1;
            e_eta01 = e_eta01 + a_beta01[mark] * (val - a_xbar01[mark]);
          }
        }
      }
      else {
        mark3 += 2;
      }

      //----- muslope and any interactions

      mark2 += 1;
      if (has_assoc01[5,m] == 1 || has_assoc01[12,m] == 1) {
        reject("muslope association structure has been removed.")
      }

      //----- muauc

      // add muauc to event submodel eta
      if (has_assoc01[6,m] == 1) { // muauc
        vector[y_qrows_for_auc01] eta_auc_tmp; // eta at all auc qpts (for submodel m)
        vector[y_qrows_for_auc01] mu_auc_tmp;  // mu  at all auc qpts (for submodel m)
        vector[y_qrows[m]] val;         // mu  following summation over auc qpts
        if (m == 1) {
          int bMat1_colshift = 0;
          int bMat2_colshift = 0;
          eta_auc_tmp = evaluate_eta01(y1_x_auc_01,
                                     y1_z1_auc_01,
                                     y1_z2_auc_01,
                                     y1_z1_id_auc_01,
                                     y1_z2_id_auc_01,
                                     yGamma1,
                                     yBeta1,
                                     bMat1,
                                     bMat2,
                                     bMat1_colshift,
                                     bMat2_colshift,
                                     intercept_type[1]);
        }
        else if (m == 2) {
          int bMat1_colshift = bK1_len[1];
          int bMat2_colshift = bK2_len[1];
          eta_auc_tmp = evaluate_eta01(y2_x_auc_01,
                                     y2_z1_auc_01,
                                     y2_z2_auc_01,
                                     y2_z1_id_auc_01,
                                     y2_z2_id_auc_01,
                                     yGamma2,
                                     yBeta2,
                                     bMat1,
                                     bMat2,
                                     bMat1_colshift,
                                     bMat2_colshift,
                                     intercept_type[2]);
        }
        else if (m == 3) {
          int bMat1_colshift = sum(bK1_len[1:2]);
          int bMat2_colshift = sum(bK2_len[1:2]);
          eta_auc_tmp = evaluate_eta01(y3_x_auc_01,
                                     y3_z1_auc_01,
                                     y3_z2_auc_01,
                                     y3_z1_id_auc_01,
                                     y3_z2_id_auc_01,
                                     yGamma3,
                                     yBeta3,
                                     bMat1,
                                     bMat2,
                                     bMat1_colshift,
                                     bMat2_colshift,
                                     intercept_type[3]);
        }
        mu_auc_tmp = evaluate_mu(eta_auc_tmp, family[m], link[m]);
        mark += 1;
        for (r in 1:y_qrows[m]) {
          vector[auc_qnodes01] val_tmp;
          vector[auc_qnodes01] wgt_tmp;
          val_tmp = mu_auc_tmp[((r-1) * auc_qnodes01 + 1):(r * auc_qnodes01)];
          wgt_tmp = auc_qwts  [((r-1) * auc_qnodes01 + 1):(r * auc_qnodes01)];
          val[r] = sum(wgt_tmp .* val_tmp);
        }
        e_eta01 += a_beta01[mark] * (val - a_xbar01[mark]);

            */
      }

    }

    //-----  shared random effects

    if (sum_size_which_b01 > 0) {
      reject("shared_b has been removed.")
    }
    if (sum_size_which_coef01 > 0) {
      reject("shared_coef has been removed.")
    }
