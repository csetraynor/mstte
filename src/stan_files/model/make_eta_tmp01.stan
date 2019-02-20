vector[y_qrows01[m]] eta_tmp01;
if (m == 1) {
  int bMat1_colshift = 0;
  int bMat2_colshift = 0;
  eta_tmp01 = evaluate_eta(y1_x_eta_01,
                         y1_z1_eta_01,
                         y1_z2_eta_01,
                         y1_z1_id_eta_01,
                         y1_z2_id_eta_01,
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
  eta_tmp01 = evaluate_eta(y2_x_eta_01,
                         y2_z1_eta_01,
                         y2_z2_eta_01,
                         y2_z1_id_eta_01,
                         y2_z2_id_eta_01,
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
  eta_tmp01 = evaluate_eta(y3_x_eta_01,
                         y3_z1_eta_01,
                         y3_z2_eta_01,
                         y3_z1_id_eta_01,
                         y3_z2_id_eta_01,
                         yGamma3,
                         yBeta3,
                         bMat1,
                         bMat2,
                         bMat1_colshift,
                         bMat2_colshift,
                         intercept_type[3]);
}
