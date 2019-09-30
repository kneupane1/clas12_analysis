/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram() {
  if (getenv("CLAS12_E") != NULL) {
    if (atof(getenv("CLAS12_E")) < 3) {
      q2_max = 1.0;
      w_max = 3.5;
      p_max = 3.0;
    } else if (atof(getenv("CLAS12_E")) < 7.9) {
      q2_max = 3.5;
      w_max = 4.0;
      p_max = 8.0;
    }
  }

  // Kinematics
  Deltat_P_after_mmsq_cut =
      new TH2D("Deltat_P_after_mmsq_cut", "Deltat_P_after_mmsq_cut", bins, p_min, p_max, bins, -10, 10);
  Deltat_pip_after_mmsq_cut =
      new TH2D("Deltat_pip_after_mmsq_cut", "Deltat_pip_after_mmsq_cut", bins, p_min, p_max, bins, -10, 10);

  Deltat_ctof_P_after_mmsq_cut =
      new TH2D("Deltat_ctof_P_after_mmsq_cut", "Deltat_ctof_P_after_mmsq_cut", bins, p_min, p_max, bins, -10, 10);
  Deltat_ctof_pip_after_mmsq_cut =
      new TH2D("Deltat_ctof_pip_after_mmsq_cut", "Deltat_ctof_pip_after_mmsq_cut", bins, p_min, p_max, bins, -10, 10);

  Deltat_P_after_mmsq_cut_singlepip =
      new TH2D("Deltat_P_after_mmsq_cut_Xpip", "Deltat_P_after_mmsq_cut_Xpip", bins, p_min, p_max, bins, -10, 10);
  Deltat_pip_after_mmsq_cut_singlepip =
      new TH2D("Deltat_pip_after_mmsq_cut_Xpip", "Deltat_pip_after_mmsq_cut_Xpip", bins, p_min, p_max, bins, -10, 10);

  Deltat_ctof_P_after_mmsq_cut_singlepip = new TH2D(
      "Deltat_ctof_P_after_mmsq_cut_Xpip", "Deltat_ctof_P_after_mmsq_cut_Xpip", bins, p_min, p_max, bins, -10, 10);
  Deltat_ctof_pip_after_mmsq_cut_singlepip = new TH2D(
      "Deltat_ctof_pip_after_mmsq_cut_Xpip", "Deltat_ctof_pip_after_mmsq_cut_Xpip", bins, p_min, p_max, bins, -10, 10);

  momentum = new TH1D("mom", "mom", bins, p_min, p_max);
  pid_size1 = new TH1D("pid_size", "pid_size", bins, -1, 10);
  pid_check = new TH1D("pid", "pid", bins, -350, 2300);
  mom_diff_e_pim = new TH1D("mom_diff_e_pim", "mom_diff_e_pim", bins, -2.2, 7.8);
  mom_diff_p_pip = new TH1D("mom_diff_p_pip", "mom_diff_p_pip", bins, -3.0, 3.0);
  mom_diff_e_pim_x = new TH1D("mom_diff_e_pim_x", "mom_diff_e_pim_x", bins, -4.0, 4.0);
  mom_diff_p_pip_x = new TH1D("mom_diff_p_pip_x", "mom_diff_p_pip_x", bins, -2.2, 2.2);
  mom_diff_e_pim_y = new TH1D("mom_diff_e_pim_y", "mom_diff_e_pim_y", bins, -4.2, 4.2);
  mom_diff_p_pip_y = new TH1D("mom_diff_p_pip_y", "mom_diff_p_pip_y", bins, -2.2, 2.2);
  mom_diff_e_pim_z = new TH1D("mom_diff_e_pim_z", "mom_diff_e_pim_z", bins, -2.2, 7.8);
  mom_diff_p_pip_z = new TH1D("mom_diff_p_pip_z", "mom_diff_p_pip_z", bins, -3.0, 3.0);
  mom_e = new TH1D("mom_e", "mom_e", bins, p_min, 7.8);
  mom_e_z = new TH1D("mom_e_z", "mom_e_z", bins, p_min, 7.8);
  mom_p = new TH1D("mom_p", "mom_p", bins, p_min, 4.2);
  mom_p_z = new TH1D("mom_p_z", "mom_p_z", bins, p_min, 4.2);
  mom_pip = new TH1D("mom_pip", "mom_pip", bins, p_min, 4.2);
  mom_pip_z = new TH1D("mom_pip_z", "mom_pip_z", bins, p_min, 4.2);
  mom_pim = new TH1D("mom_pim", "mom_pim", bins, p_min, 4.2);
  mom_pim_z = new TH1D("mom_pim_z", "mom_pim_z", bins, -1, 3.2);

  vertex_vz = new TH1D("vertex_position", "vertex_position", bins, -40, 40);
  W_thrown = new TH1D("W_thrown", "W_thrown", bins, p_min, 5);
  W_vs_Q2_thrown = new TH2D("W_vs_Q2_thrown", "W_vs_Q2_thrown", bins, 1, 3.5, bins, 0.7, 3.6);
  W_vs_Q2_all_sec = new TH2D("W_vs_Q2_all_sec", "W_vs_Q2_all_sec", bins, p_min, 5, bins, p_min, 8.5);
  lu_side_distribution = new TH1D("lu_side_distribution", "lu_side_distribution", 50, 0, 400);
  lv_side_distribution = new TH1D("lv_side_distribution", "lv_side_distribution", 50, 0, 450);
  lw_side_distribution = new TH1D("lw_side_distribution", "lw_side_distribution", 50, 0, 450);

  W_hist_Xpip_all_sec = new TH1D("W_Xpip_all_sec", "W_Xpip_all_sec", bins, 0.9, 2.8);

  // makeHists_deltat();
  // makeHists_MomVsBeta();
  makeHists_WvsQ2();
  makeHists_Angles();
  makeHists_MM();
  // Make_hist_cc();
  makeHists_EC_sf();
  makeHists_sf_slices();
  makeHists_pcal_fid_cuts();
}

Histogram::~Histogram() {}
// W and Q^2

float Histogram::mm_lim_min(int mm_number, int mm_events_number) {
  if (mm_number == 0 && mm_events_number < 4) {
    return -2.25;
  } else if (mm_number == 0 && mm_events_number > 4) {
    return -1.0;
  } else if (mm_number == 0 && mm_events_number == 4) {
    return 0;
  } else if (mm_number == 1 && mm_events_number < 4 && mm_events_number != 1) {
    return -1.5;  // 0.05 original value
  } else if (mm_number == 1 && mm_events_number > 4) {
    return -2;
  } else if (mm_number == 1 && mm_events_number == 1) {
    return -0.20;
  } else if (mm_number == 1 && mm_events_number == 4) {
    return 0;
  } else {
    return -5;
  }
}
float Histogram::mm_lim_max(int mm_number, int mm_events_number) {
  if (mm_number == 0 && mm_events_number < 4) {
    return 3.35;
  } else if (mm_number == 0 && mm_events_number > 4) {
    return 4.50;
  } else if (mm_number == 0 && mm_events_number == 4) {
    return 6.0;
  } else if (mm_number == 1 && mm_events_number < 4 && mm_events_number != 1) {
    return 10.5;  // 0.05 original value
  } else if (mm_number == 1 && mm_events_number == 1) {
    return 0.20;
  } else if (mm_number == 1 && mm_events_number > 4) {
    return 19.0;
  } else if (mm_number == 1 && mm_events_number == 4) {
    return 20;
  } else {
    return 20;
  }
}
float Histogram::w_hist_lim_max(int w_num) {
  if (w_num == 0) {
    return 2.7;
  } else if (w_num == 1) {
    return 3.6;
  } else {
    return 3.6;
  }
}

void Histogram::makeHists_pcal_fid_cuts() {
  for (int j = 0; j < cut_y_n; j++) {
    hname.append("PCAL_FID_CUTS");
    htitle.append("PCAL_FID_CUTS");
    hname.append("_");
    htitle.append("_");
    hname.append(cut_without_cut_name[j]);
    htitle.append(cut_without_cut_name[j]);
    PCAL_FID_CUT[j] = new TH2D(hname.c_str(), htitle.c_str(), bins, -400, 400, bins, -400, 400);
    hname.clear();
    htitle.clear();

    hname.append("DCr1_FID_CUTS");
    htitle.append("DCr1_FID_CUTS");
    hname.append("_");
    htitle.append("_");
    hname.append(cut_without_cut_name[j]);
    htitle.append(cut_without_cut_name[j]);
    DCr1_FID_CUT[j] = new TH2D(hname.c_str(), htitle.c_str(), bins, -155, 155, bins, -155, 155);
    hname.clear();
    htitle.clear();

    hname.append("DCr2_FID_CUTS");
    htitle.append("DCr2_FID_CUTS");
    hname.append("_");
    htitle.append("_");
    hname.append(cut_without_cut_name[j]);
    htitle.append(cut_without_cut_name[j]);
    DCr2_FID_CUT[j] = new TH2D(hname.c_str(), htitle.c_str(), bins, -200, 200, bins, -200, 200);
    hname.clear();
    htitle.clear();

    hname.append("DCr3_FID_CUTS");
    htitle.append("DCr3_FID_CUTS");
    hname.append("_");
    htitle.append("_");
    hname.append(cut_without_cut_name[j]);
    htitle.append(cut_without_cut_name[j]);
    DCr3_FID_CUT[j] = new TH2D(hname.c_str(), htitle.c_str(), bins, -300, 300, bins, -300, 300);
    hname.clear();
    htitle.clear();
  }
}

void Histogram::Fill_hist_PCAL_FID_CUT(float x_PCAL, float y_PCAL) { PCAL_FID_CUT[0]->Fill(x_PCAL, y_PCAL); }
void Histogram::Fill_hist_PCAL_without_FID_CUT(float x_PCAL, float y_PCAL) { PCAL_FID_CUT[1]->Fill(x_PCAL, y_PCAL); }
void Histogram::Fill_hist_DC_FID_CUT(float R1X, float R1Y, float R2X, float R2Y, float R3X, float R3Y) {
  DCr1_FID_CUT[0]->Fill(R1X, R1Y);
  DCr2_FID_CUT[0]->Fill(R2X, R2Y);
  DCr3_FID_CUT[0]->Fill(R3X, R3Y);
}
void Histogram::Fill_hist_DC_without_FID_CUT(float R1X, float R1Y, float R2X, float R2Y, float R3X, float R3Y) {
  DCr1_FID_CUT[1]->Fill(R1X, R1Y);
  DCr2_FID_CUT[1]->Fill(R2X, R2Y);
  DCr3_FID_CUT[1]->Fill(R3X, R3Y);
}

void Histogram::makeHists_EC_sf() {
  for (int i = 0; i < sec_num; i++) {
    hname.append("EC_sampling_fraction");
    htitle.append("EC_sampling_fraction");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    EC_sampling_fraction[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, zero, 0.4);
    hname.clear();
    htitle.clear();

    hname.append("PCAL_VS_ECAL");
    htitle.append("PCAL_VS_ECAL");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    PCAL_VS_ECAL[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, zero, 1.5, bins, zero, 1.5);
    hname.clear();
    htitle.clear();
  }
}

void Histogram::makeHists_Angles() {
  for (int i = 0; i < w_range_num; i++) {
    for (int j = 0; j < q2_range_num; j++) {
      hname.append("theta_P_thrown");
      htitle.append("theta_P_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_P_cm_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("theta_pip_thrown");
      htitle.append("theta_pip_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pip_cm_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("theta_pim_thrown");
      htitle.append("theta_pim_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pim_cm_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("phi_P_thrown");
      htitle.append("phi_P_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      phi_P_cm_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("phi_pip_thrown");
      htitle.append("phi_pip_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      phi_pip_cm_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("phi_pim_thrown");
      htitle.append("phi_pim_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      phi_pim_cm_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("alpha_P_cm_thrown");
      htitle.append("alpha_P_cm_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      alpha_P_cm_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 360);
      hname.clear();
      htitle.clear();
      hname.append("alpha_pip_cm_thrown");
      htitle.append("alpha_pip_cm_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      alpha_pip_cm_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 360);
      hname.clear();
      htitle.clear();
      hname.append("alpha_pim_cm_thrown");
      htitle.append("alpha_pim_cm_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      alpha_pim_cm_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 360);
      hname.clear();
      htitle.clear();

      hname.append("theta_P");
      htitle.append("theta_P");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_P_cm[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("theta_pip");
      htitle.append("theta_pip");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pip_cm[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("theta_pim");
      htitle.append("theta_pim");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pim_cm[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("phi_P");
      htitle.append("phi_P");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      phi_P_cm[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("phi_pip");
      htitle.append("phi_pip");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      phi_pip_cm[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("phi_pim");
      htitle.append("phi_pim");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      phi_pim_cm[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("alpha_P_cm");
      htitle.append("alpha_P_cm");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      alpha_P_cm[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 360);
      hname.clear();
      htitle.clear();
      hname.append("alpha_pip_cm");
      htitle.append("alpha_pip_cm");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      alpha_pip_cm[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 360);
      hname.clear();
      htitle.clear();
      hname.append("alpha_pim_cm");
      htitle.append("alpha_pim_cm");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      alpha_pim_cm[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 360);
      hname.clear();
      htitle.clear();
    }
  }
}
void Histogram::makeHists_WvsQ2() {
  for (int i = 0; i < w_range_num; i++) {
    for (int j = 0; j < q2_range_num; j++) {
      hname.append("theta_p_vs_inv_mass_pip_pim_thrown");
      htitle.append("theta_p_vs_inv_mass_pip_pim_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);

      theta_P_vs_mass_pip_pim_thrown[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 0.2, Histogram::w_hist_lim_max(i), 100, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("theta_pim_vs_inv_mass_Ppip_thrown");
      htitle.append("theta_pim_vs_inv_mass_Ppip_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pim_vs_mass_Ppip_thrown[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i), 100, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("theta_pip_vs_inv_mass_Ppim_thrown");
      htitle.append("theta_pip_vs_inv_mass_Ppim_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pip_vs_mass_Ppim_thrown[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i), 100, 0, 180);
      hname.clear();
      htitle.clear();

      hname.append("theta_P_lab_vs_inv_mass_pip_pim_thrown");
      htitle.append("theta_P_lab_vs_inv_mass_pip_pim_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_P_lab_vs_mass_pip_pim_thrown[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 0.2, Histogram::w_hist_lim_max(i), 100, 0, 75);
      hname.clear();
      htitle.clear();

      hname.append("theta_pim_lab_vs_inv_mass_P_pip_thrown");
      htitle.append("theta_pim_lab_vs_inv_mass_P_pip_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pim_lab_vs_mass_Ppip_thrown[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i), 100, 0, 125);
      hname.clear();
      htitle.clear();

      hname.append("theta_pip_lab_vs_inv_mass_P_pim_thrown");
      htitle.append("theta_pip_lab_vs_inv_mass_P_pim_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pip_lab_vs_mass_Ppim_thrown[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i), 100, 0, 125);
      hname.clear();
      htitle.clear();

      hname.append("W_P#pi+#pi-_thrown");
      htitle.append("W_P#pi+#pi-_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      W_hist_twopi_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i));
      hname.clear();
      htitle.clear();

      hname.append("invariant_mass P#pi+_thrown");
      htitle.append("invariant_mass P#pi+_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      inv_mass_P_pip_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i));
      hname.clear();
      htitle.clear();

      hname.append("invariant_mass P#pi-_thrown");
      htitle.append("invariant_mass P#pi-_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      inv_mass_P_pim_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i));
      hname.clear();
      htitle.clear();

      hname.append("invariant_mass #pi+#pi-_thrown");
      htitle.append("invariant_mass #pi+#pi-_thrown");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      inv_mass_pip_pim_thrown[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0.2, Histogram::w_hist_lim_max(i));
      hname.clear();
      htitle.clear();

      hname.append("theta_p_vs_inv_mass_pip_pim");
      htitle.append("theta_p_vs_inv_mass_pip_pim");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_P_vs_mass_pip_pim[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 0.2, Histogram::w_hist_lim_max(i), 100, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("theta_pim_vs_inv_mass_Ppip");
      htitle.append("theta_pim_vs_inv_mass_Ppip");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pim_vs_mass_Ppip[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i), 100, 0, 180);
      hname.clear();
      htitle.clear();
      hname.append("theta_pip_vs_inv_mass_Ppim");
      htitle.append("theta_pip_vs_inv_mass_Ppim");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append(" ");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pip_vs_mass_Ppim[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i), 100, 0, 180);
      hname.clear();
      htitle.clear();

      hname.append("theta_P_lab_vs_inv_mass_pip_pim");
      htitle.append("theta_P_lab_vs_inv_mass_pip_pim");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_P_lab_vs_mass_pip_pim[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 0.2, Histogram::w_hist_lim_max(i), 100, 0, 75);
      hname.clear();
      htitle.clear();

      hname.append("theta_pim_lab_vs_inv_mass_P_pip");
      htitle.append("theta_pim_lab_vs_inv_mass_P_pip");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pim_lab_vs_mass_Ppip[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i), 100, 0, 125);
      hname.clear();
      htitle.clear();

      hname.append("theta_pip_lab_vs_inv_mass_P_pim");
      htitle.append("theta_pip_lab_vs_inv_mass_P_pim");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      theta_pip_lab_vs_mass_Ppim[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i), 100, 0, 125);
      hname.clear();
      htitle.clear();

      // hname.append("W_vs_Q2_in_range");
      // htitle.append("W_vs_Q2_in_range");
      // hname.append("_W_");
      // htitle.append(" W ");
      // hname.append(w_range_name[i]);
      // htitle.append(w_range_name[i]);
      // htitle.append("_Q2_");
      // htitle.append(" Q2 ");
      // hname.append(q2_range_name[j]);
      // htitle.append(q2_range_name[j]);
      // W_vs_Q2_all_sec[i][j] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, 5, bins, p_min, q2_max);
      // hname.clear();
      // htitle.clear();

      hname.append("Prot_mass_w_vs_Q2");
      htitle.append("Prot_mass_w_vs_Q2");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      Prot_mass_w_vs_Q2[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, Histogram::w_hist_lim_max(i), bins, p_min, q2_max);
      hname.clear();
      htitle.clear();

      hname.append("Prot_theta_lab_vs_W");
      htitle.append("Prot_theta_lab_vs_W");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      Prot_theta_lab_vs_W[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, Histogram::w_hist_lim_max(i), 100, 0, 75);
      hname.clear();
      htitle.clear();

      hname.append("Pip_mass_w_vs_Q2");
      htitle.append("Pip_mass_w_vs_Q2");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      Pip_mass_w_vs_Q2[i][j] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, 1, bins, p_min, q2_max);
      hname.clear();
      htitle.clear();

      hname.append("theta_Pip_lab_vs_w");
      htitle.append("theta_Pip_lab_vs_w");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      Pip_theta_lab_vs_W[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, Histogram::w_hist_lim_max(i), 100, 0, 120);
      hname.clear();
      htitle.clear();

      hname.append("Pim_theta_lab_vs_W");
      htitle.append("Pim_theta_lab_vs_W");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      Pim_theta_lab_vs_W[i][j] =
          new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, Histogram::w_hist_lim_max(i), 100, 0, 120);
      hname.clear();
      htitle.clear();

      hname.append("W_P#pi+#pi-");
      htitle.append("W_P#pi+#pi-");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      W_hist_twopi_all_sec[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i));
      hname.clear();
      htitle.clear();

      hname.append("invariant_mass P#pi+");
      htitle.append("invariant_mass P#pi+");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      inv_mass_P_pip_all_sec[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i));
      hname.clear();
      htitle.clear();

      hname.append("invariant_mass P#pi-");
      htitle.append("invariant_mass P#pi-");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      inv_mass_P_pim_all_sec[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 1, Histogram::w_hist_lim_max(i));
      hname.clear();
      htitle.clear();

      hname.append("invariant_mass #pi+#pi-");
      htitle.append("invariant_mass #pi+#pi-");
      hname.append("_");
      htitle.append(" ");
      hname.append(w_range_name[i]);
      htitle.append(w_range_name[i]);
      htitle.append("_");
      htitle.append(" ");
      hname.append(q2_range_name[j]);
      htitle.append(q2_range_name[j]);
      inv_mass_pip_pim_all_sec[i][j] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0.2, Histogram::w_hist_lim_max(i));
      hname.clear();
      htitle.clear();
    }
  }
  for (int i = 0; i < sec_num; i++) {
    hname.append("WvsQ2");
    htitle.append("WvsQ2");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    W_vs_Q2[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, w_max, bins, p_min, q2_max);
    hname.clear();
    htitle.clear();

    hname.append("W_hist");
    htitle.append("W_hist");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    W_hist[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, p_min, w_max);
    hname.clear();
    htitle.clear();

    hname.append("invariant_mass_ep");
    htitle.append("invariant_mass_ep");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    W_hist_ep[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, p_min, w_max);
    hname.clear();
    htitle.clear();

    hname.append("W_P#pi+#pi-");
    htitle.append("W_P#pi+#pi-");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    W_hist_2pi[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 1, 3.8);
    hname.clear();
    htitle.clear();

    hname.append("invariant_mass P#pi+");
    htitle.append("invariant_mass P#pi+");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    W_hist_delta_pp[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 1, 3.1);
    hname.clear();
    htitle.clear();

    hname.append("invariant_mass P#pi-");
    htitle.append("invariant_mass P#pi-");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    W_hist_delta_zero[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 1, 3.1);
    hname.clear();
    htitle.clear();

    hname.append("invariant_mass #pi+#pi-");
    htitle.append("invariant_mass #pi+#pi-");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    W_hist_rho[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0.2, 2);
    hname.clear();
    htitle.clear();

    hname.append("invariant_mass N#pi+ ");
    htitle.append("invariant_mass N#pi+");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    W_hist_singlepip[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0.7, 3.1);
    hname.clear();
    htitle.clear();

    hname.append("Q2_hist");
    htitle.append("Q2_hist");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    Q2_hist[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, p_min, q2_max);
    hname.clear();
    htitle.clear();
    for (int j = 0; j < cut_y_n; j++) {
      hname.append("Wvs_mmSQ_e(p,p'X)e'");
      htitle.append("Wvs_mmSQ_e(p,p'X)e'");
      hname.append("_");
      htitle.append(" ");
      hname.append(sec_name[i]);
      htitle.append(sec_name[i]);
      hname.append("_");
      htitle.append(" ");
      hname.append(cut_name[j]);
      htitle.append(cut_name[j]);
      W_vs_mmSQ_ep[i][j] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, w_max, bins, -4, 4);
      hname.clear();
      htitle.clear();
      hname.append("Wvs_mmSQ_e(p,p'pi+pi-X)e'");
      htitle.append("Wvs_mmSQ_e(p,p'pi+pi-X)e'");
      hname.append("_");
      htitle.append(" ");
      hname.append(sec_name[i]);
      htitle.append(sec_name[i]);
      hname.append("_");
      htitle.append(" ");
      hname.append(cut_name[j]);
      htitle.append(cut_name[j]);
      W_vs_mmSQ_2pi[i][j] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, w_max, bins, -5, 5);
      hname.clear();
      htitle.clear();
      hname.clear();
      htitle.clear();
      hname.append("Wvs_mmSQ_e(p,pi+X)e'");
      htitle.append("Wvs_mmSQ_e(p,pi+X)e'");
      hname.append("_");
      htitle.append(" ");
      hname.append(sec_name[i]);
      htitle.append(sec_name[i]);
      hname.append("_");
      htitle.append(" ");
      hname.append(cut_name[j]);
      htitle.append(cut_name[j]);
      W_vs_mmSQ_singlepip[i][j] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, w_max, bins, -5, 5);
      hname.clear();
      htitle.clear();
    }
  }
}

// void Histogram::Fill_WvsQ2_range(float W, float Q2, float weight) {
//   if (W < 2.0)
//     W_vs_Q2[0]->Fill(W, Q2, weight);
//   else if (W > 2.0 && W < 2.5)
//     W_vs_Q2[1]->Fill(W, Q2, weight);
//   else if (W > 2.5 && W < 3.0)
//     W_vs_Q2[2]->Fill(W, Q2, weight);
//   else if (W > 3.0 && W < 3.5)
//     W_vs_Q2[3]->Fill(W, Q2, weight);
// }
void Histogram::Fill_WvsQ2(float W, float Q2, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_vs_Q2[sec_number]->Fill(W, Q2, weight);
    W_hist[sec_number]->Fill(W, weight);

    Q2_hist[sec_number]->Fill(Q2, weight);
  }
}
void Histogram::Fill_EC_sampling_fraction(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    EC_sampling_fraction[sec_number]->Fill(momentum, sf, weight);
  }
}

void Histogram::makeHists_sf_slices() {
  for (int i = 0; i < sec_num; i++) {
    for (int n = 0; n < slice_num; n++) {
      hname.append("sampling_fraction_slice");
      htitle.append("sampling_fraction_slice");
      hname.append("_");
      htitle.append(" ");
      hname.append(sec_name[i]);
      htitle.append(sec_name[i]);
      hname.append("_");
      htitle.append(" ");
      hname.append(slice_name[n]);
      htitle.append(slice_name[n]);
      sf_[i][n] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 0.5);
      hname.clear();
      htitle.clear();
    }
  }
}
void Histogram::Fill_1d_sampling_fraction_1(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    // sf_[sec_number][1] -> Fit("gaus", "", "", 0.18, 0.28);
    sf_[sec_number][1]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_2(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][2]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_3(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][3]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_4(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][4]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_5(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][5]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_6(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][6]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_7(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][7]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_8(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][8]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_9(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][9]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_10(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][10]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_11(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][11]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_12(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][12]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_13(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][13]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_14(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][14]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_15(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][15]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_16(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][16]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_17(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][17]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_18(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][18]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_19(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][19]->Fill(sf, weight);
  }
}
void Histogram::Fill_1d_sampling_fraction_20(float momentum, float sf, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    sf_[sec_number][20]->Fill(sf, weight);
  }
}

void Histogram::Fill_PCAL_VS_ECAL(float pcal, float ecal, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    PCAL_VS_ECAL[sec_number]->Fill(pcal, ecal, weight);
  }
}
void Histogram::Fill_WvsmmSQ_ep(float W, float mmSQ, int sec_number, float weight) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_hist_ep[sec_number]->Fill(W, weight);
    W_vs_mmSQ_ep[sec_number][0]->Fill(W, mmSQ, weight);
  }
}
void Histogram::Fill_W_2pi_all_sec(float W, float q2, float W_2pi, float W_dpp, float delta_zero_, float rho_,
                                   float weight) {
  if (W < 2.5) {
    if (q2 < 4.5) {
      W_hist_twopi_all_sec[0][0]->Fill(W_2pi, weight);
      inv_mass_P_pip_all_sec[0][0]->Fill(W_dpp, weight);
      inv_mass_P_pim_all_sec[0][0]->Fill(delta_zero_, weight);
      inv_mass_pip_pim_all_sec[0][0]->Fill(rho_, weight);
    } else if (q2 > 4.5 && q2 < 15.0) {
      W_hist_twopi_all_sec[0][1]->Fill(W_2pi, weight);
      inv_mass_P_pip_all_sec[0][1]->Fill(W_dpp, weight);
      inv_mass_P_pim_all_sec[0][1]->Fill(delta_zero_, weight);
      inv_mass_pip_pim_all_sec[0][1]->Fill(rho_, weight);
      // } else if (q2 > 5.0 && q2 < 15.0) {
      //   W_hist_twopi_all_sec[0][2]->Fill(W_2pi, weight);
      //   inv_mass_P_pip_all_sec[0][2]->Fill(W_dpp, weight);
      //   inv_mass_P_pim_all_sec[0][2]->Fill(delta_zero_, weight);
      //   inv_mass_pip_pim_all_sec[0][2]->Fill(rho_, weight);
    }
  } else if (W > 2.5 && W < 3.5) {
    if (q2 < 4.5) {
      W_hist_twopi_all_sec[1][0]->Fill(W_2pi, weight);
      inv_mass_P_pip_all_sec[1][0]->Fill(W_dpp, weight);
      inv_mass_P_pim_all_sec[1][0]->Fill(delta_zero_, weight);
      inv_mass_pip_pim_all_sec[1][0]->Fill(rho_, weight);
    } else if (q2 > 4.5 && q2 < 15.0) {
      W_hist_twopi_all_sec[1][1]->Fill(W_2pi, weight);
      inv_mass_P_pip_all_sec[1][1]->Fill(W_dpp, weight);
      inv_mass_P_pim_all_sec[1][1]->Fill(delta_zero_, weight);
      inv_mass_pip_pim_all_sec[1][1]->Fill(rho_, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     W_hist_twopi_all_sec[1][2]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_all_sec[1][2]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_all_sec[1][2]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_all_sec[1][2]->Fill(rho_, weight);
      //   }
      // } else if (W > 2.5 && W < 3.0) {
      //   if (q2 < 2.0) {
      //     W_hist_twopi_all_sec[2][0]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_all_sec[2][0]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_all_sec[2][0]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_all_sec[2][0]->Fill(rho_, weight);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     W_hist_twopi_all_sec[2][1]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_all_sec[2][1]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_all_sec[2][1]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_all_sec[2][1]->Fill(rho_, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     W_hist_twopi_all_sec[2][2]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_all_sec[2][2]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_all_sec[2][2]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_all_sec[2][2]->Fill(rho_, weight);
      //   }
      // } else if (W > 3.0 && W < 3.5) {
      //   if (q2 < 2.0) {
      //     W_hist_twopi_all_sec[3][0]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_all_sec[3][0]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_all_sec[3][0]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_all_sec[3][0]->Fill(rho_, weight);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     W_hist_twopi_all_sec[3][1]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_all_sec[3][1]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_all_sec[3][1]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_all_sec[3][1]->Fill(rho_, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     W_hist_twopi_all_sec[3][2]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_all_sec[3][2]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_all_sec[3][2]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_all_sec[3][2]->Fill(rho_, weight);
    }
  }
}
void Histogram::Fill_W_2pi_thrown(float W, float q2, float W_2pi, float W_dpp, float delta_zero_, float rho_,
                                  float weight) {
  if (W < 2.5) {
    if (q2 < 4.5) {
      W_hist_twopi_thrown[0][0]->Fill(W_2pi, weight);
      inv_mass_P_pip_thrown[0][0]->Fill(W_dpp, weight);
      inv_mass_P_pim_thrown[0][0]->Fill(delta_zero_, weight);
      inv_mass_pip_pim_thrown[0][0]->Fill(rho_, weight);
    } else if (q2 > 4.5 && q2 < 15.0) {
      W_hist_twopi_thrown[0][1]->Fill(W_2pi, weight);
      inv_mass_P_pip_thrown[0][1]->Fill(W_dpp, weight);
      inv_mass_P_pim_thrown[0][1]->Fill(delta_zero_, weight);
      inv_mass_pip_pim_thrown[0][1]->Fill(rho_, weight);
      // } else if (q2 > 5.0 && q2 < 15.0) {
      //   W_hist_twopi_thrown[0][2]->Fill(W_2pi, weight);
      //   inv_mass_P_pip_thrown[0][2]->Fill(W_dpp, weight);
      //   inv_mass_P_pim_thrown[0][2]->Fill(delta_zero_, weight);
      //   inv_mass_pip_pim_thrown[0][2]->Fill(rho_, weight);
    }
  } else if (W > 2.5 && W < 3.5) {
    if (q2 < 4.5) {
      W_hist_twopi_thrown[1][0]->Fill(W_2pi, weight);
      inv_mass_P_pip_thrown[1][0]->Fill(W_dpp, weight);
      inv_mass_P_pim_thrown[1][0]->Fill(delta_zero_, weight);
      inv_mass_pip_pim_thrown[1][0]->Fill(rho_, weight);
    } else if (q2 > 4.5 && q2 < 15.0) {
      W_hist_twopi_thrown[1][1]->Fill(W_2pi, weight);
      inv_mass_P_pip_thrown[1][1]->Fill(W_dpp, weight);
      inv_mass_P_pim_thrown[1][1]->Fill(delta_zero_, weight);
      inv_mass_pip_pim_thrown[1][1]->Fill(rho_, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     W_hist_twopi_thrown[1][2]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_thrown[1][2]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_thrown[1][2]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_thrown[1][2]->Fill(rho_, weight);
      //   }
      // } else if (W > 2.5 && W < 3.0) {
      //   if (q2 < 2.0) {
      //     W_hist_twopi_thrown[2][0]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_thrown[2][0]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_thrown[2][0]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_thrown[2][0]->Fill(rho_, weight);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     W_hist_twopi_thrown[2][1]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_thrown[2][1]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_thrown[2][1]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_thrown[2][1]->Fill(rho_, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     W_hist_twopi_thrown[2][2]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_thrown[2][2]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_thrown[2][2]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_thrown[2][2]->Fill(rho_, weight);
      //   }
      // } else if (W > 3.0 && W < 3.5) {
      //   if (q2 < 2.0) {
      //     W_hist_twopi_thrown[3][0]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_thrown[3][0]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_thrown[3][0]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_thrown[3][0]->Fill(rho_, weight);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     W_hist_twopi_thrown[3][1]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_thrown[3][1]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_thrown[3][1]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_thrown[3][1]->Fill(rho_, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     W_hist_twopi_thrown[3][2]->Fill(W_2pi, weight);
      //     inv_mass_P_pip_thrown[3][2]->Fill(W_dpp, weight);
      //     inv_mass_P_pim_thrown[3][2]->Fill(delta_zero_, weight);
      //     inv_mass_pip_pim_thrown[3][2]->Fill(rho_, weight);
    }
  }
}

void Histogram::Fill_theta_P_inv_mass(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_P_vs_mass_pip_pim[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_vs_mass_pip_pim[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_P_vs_mass_pip_pim[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_P_vs_mass_pip_pim[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_vs_mass_pip_pim[1][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_P_vs_mass_pip_pim[1][2]->Fill(inv_mass, theta, wt);
      // }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_P_vs_mass_pip_pim[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_vs_mass_pip_pim[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_vs_mass_pip_pim[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_P_vs_mass_pip_pim[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_vs_mass_pip_pim[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_vs_mass_pip_pim[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}
void Histogram::Fill_theta_pim_inv_mass(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_pim_vs_mass_Ppip[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pim_vs_mass_Ppip[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_pim_vs_mass_Ppip[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_pim_vs_mass_Ppip[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pim_vs_mass_Ppip[1][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_vs_mass_Ppip[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_pim_vs_mass_Ppip[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pim_vs_mass_Ppip[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_vs_mass_Ppip[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_pim_vs_mass_Ppip[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pim_vs_mass_Ppip[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_vs_mass_Ppip[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}
void Histogram::Fill_theta_pip_inv_mass(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_pip_vs_mass_Ppim[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pip_vs_mass_Ppim[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_pip_vs_mass_Ppim[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_pip_vs_mass_Ppim[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.50 && q2 < 15.0) {
      theta_pip_vs_mass_Ppim[1][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_vs_mass_Ppim[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_pip_vs_mass_Ppim[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pip_vs_mass_Ppim[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_vs_mass_Ppim[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_pip_vs_mass_Ppim[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pip_vs_mass_Ppim[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_vs_mass_Ppim[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}
void Histogram::Fill_theta_P_lab_inv_mass(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_P_lab_vs_mass_pip_pim[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_lab_vs_mass_pip_pim[0][1]->Fill(inv_mass, theta, wt);
    } else if (q2 > 5.0 && q2 < 10.0) {
      theta_P_lab_vs_mass_pip_pim[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 5.5) {
    if (q2 < 4.5) {
      theta_P_lab_vs_mass_pip_pim[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_lab_vs_mass_pip_pim[1][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_lab_vs_mass_pip_pim[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_P_lab_vs_mass_pip_pim[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_lab_vs_mass_pip_pim[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_lab_vs_mass_pip_pim[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_P_lab_vs_mass_pip_pim[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_lab_vs_mass_pip_pim[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_lab_vs_mass_pip_pim[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}
void Histogram::Fill_theta_pim_lab_inv_mass(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_pim_lab_vs_mass_Ppip[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pim_lab_vs_mass_Ppip[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_pim_lab_vs_mass_Ppip[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_pim_lab_vs_mass_Ppip[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      //     theta_pim_lab_vs_mass_Ppip[1][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_lab_vs_mass_Ppip[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_pim_lab_vs_mass_Ppip[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pim_lab_vs_mass_Ppip[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_lab_vs_mass_Ppip[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_pim_lab_vs_mass_Ppip[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pim_lab_vs_mass_Ppip[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_lab_vs_mass_Ppip[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}
void Histogram::Fill_theta_pip_lab_inv_mass(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_pip_lab_vs_mass_Ppim[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pip_lab_vs_mass_Ppim[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_pip_lab_vs_mass_Ppim[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_pip_lab_vs_mass_Ppim[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pip_lab_vs_mass_Ppim[1][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_lab_vs_mass_Ppim[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_pip_lab_vs_mass_Ppim[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pip_lab_vs_mass_Ppim[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_lab_vs_mass_Ppim[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_pip_lab_vs_mass_Ppim[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pip_lab_vs_mass_Ppim[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_lab_vs_mass_Ppim[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}

void Histogram::Fill_theta_P_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_P_vs_mass_pip_pim_thrown[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_vs_mass_pip_pim_thrown[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_P_vs_mass_pip_pim_thrown[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_P_vs_mass_pip_pim_thrown[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_vs_mass_pip_pim_thrown[1][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_vs_mass_pip_pim_thrown[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_P_vs_mass_pip_pim_thrown[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_vs_mass_pip_pim_thrown[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_vs_mass_pip_pim_thrown[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_P_vs_mass_pip_pim_thrown[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_vs_mass_pip_pim_thrown[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_vs_mass_pip_pim_thrown[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}
void Histogram::Fill_theta_pim_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_pim_vs_mass_Ppip_thrown[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pim_vs_mass_Ppip_thrown[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_pim_vs_mass_Ppip_thrown[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_pim_vs_mass_Ppip_thrown[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pim_vs_mass_Ppip_thrown[1][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_pim_vs_mass_Ppip_thrown[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_pim_vs_mass_Ppip_thrown[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pim_vs_mass_Ppip_thrown[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_vs_mass_Ppip_thrown[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_pim_vs_mass_Ppip_thrown[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pim_vs_mass_Ppip_thrown[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_vs_mass_Ppip_thrown[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}
void Histogram::Fill_theta_pip_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_pip_vs_mass_Ppim_thrown[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pip_vs_mass_Ppim_thrown[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_pip_vs_mass_Ppim_thrown[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_pip_vs_mass_Ppim_thrown[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pip_vs_mass_Ppim_thrown[1][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_vs_mass_Ppim_thrown[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_pip_vs_mass_Ppim_thrown[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pip_vs_mass_Ppim_thrown[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_vs_mass_Ppim_thrown[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_pip_vs_mass_Ppim_thrown[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pip_vs_mass_Ppim_thrown[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_vs_mass_Ppim_thrown[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}
void Histogram::Fill_theta_P_lab_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_P_lab_vs_mass_pip_pim_thrown[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_lab_vs_mass_pip_pim_thrown[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_P_lab_vs_mass_pip_pim_thrown[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_P_lab_vs_mass_pip_pim_thrown[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_lab_vs_mass_pip_pim_thrown[1][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_lab_vs_mass_pip_pim_thrown[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_P_lab_vs_mass_pip_pim_thrown[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_lab_vs_mass_pip_pim_thrown[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_lab_vs_mass_pip_pim_thrown[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_P_lab_vs_mass_pip_pim_thrown[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_lab_vs_mass_pip_pim_thrown[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_P_lab_vs_mass_pip_pim_thrown[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}
void Histogram::Fill_theta_pim_lab_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_pim_lab_vs_mass_Ppip_thrown[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pim_lab_vs_mass_Ppip_thrown[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_pim_lab_vs_mass_Ppip_thrown[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_pim_lab_vs_mass_Ppip_thrown[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pim_lab_vs_mass_Ppip_thrown[1][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_lab_vs_mass_Ppip_thrown[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_pim_lab_vs_mass_Ppip_thrown[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pim_lab_vs_mass_Ppip_thrown[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_lab_vs_mass_Ppip_thrown[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_pim_lab_vs_mass_Ppip_thrown[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pim_lab_vs_mass_Ppip_thrown[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pim_lab_vs_mass_Ppip_thrown[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}
void Histogram::Fill_theta_pip_lab_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_pip_lab_vs_mass_Ppim_thrown[0][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pip_lab_vs_mass_Ppim_thrown[0][1]->Fill(inv_mass, theta, wt);
      // } else if (q2 > 5.0 && q2 < 10.0) {
      //   theta_pip_lab_vs_mass_Ppim_thrown[0][2]->Fill(inv_mass, theta, wt);
    }
  } else if (W >= 2.5 && W <= 3.5) {
    if (q2 < 4.5) {
      theta_pip_lab_vs_mass_Ppim_thrown[1][0]->Fill(inv_mass, theta, wt);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_pip_lab_vs_mass_Ppim_thrown[1][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_lab_vs_mass_Ppim_thrown[1][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 2.5 && W <= 3.0) {
      //   if (q2 < 2.0) {
      //     theta_pip_lab_vs_mass_Ppim_thrown[2][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pip_lab_vs_mass_Ppim_thrown[2][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_lab_vs_mass_Ppim_thrown[2][2]->Fill(inv_mass, theta, wt);
      //   }
      // } else if (W >= 3.0 && W <= 3.5) {
      //   if (q2 < 2.0) {
      //     theta_pip_lab_vs_mass_Ppim_thrown[3][0]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_pip_lab_vs_mass_Ppim_thrown[3][1]->Fill(inv_mass, theta, wt);
      //   } else if (q2 > 5.0 && q2 < 10.0) {
      //     theta_pip_lab_vs_mass_Ppim_thrown[3][2]->Fill(inv_mass, theta, wt);
    }
  }
}

void Histogram::Fill_Angles(float W, float q2, float theta_p, float theta_pip, float theta_pim, float phi_p,
                            float phi_pip, float phi_pim, float alpha_p, float alpha_pip, float alpha_pim,
                            float weight) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_P_cm[0][0]->Fill(theta_p, weight);
      theta_pip_cm[0][0]->Fill(theta_pip, weight);
      theta_pim_cm[0][0]->Fill(theta_pim, weight);

      phi_P_cm[0][0]->Fill(phi_p, weight);
      phi_pip_cm[0][0]->Fill(phi_pip, weight);
      phi_pim_cm[0][0]->Fill(phi_pim, weight);

      alpha_P_cm[0][0]->Fill(alpha_p, weight);
      alpha_pip_cm[0][0]->Fill(alpha_pip, weight);
      alpha_pim_cm[0][0]->Fill(alpha_pim, weight);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_cm[0][1]->Fill(theta_p, weight);
      theta_pip_cm[0][1]->Fill(theta_pip, weight);
      theta_pim_cm[0][1]->Fill(theta_pim, weight);

      phi_P_cm[0][1]->Fill(phi_p, weight);
      phi_pip_cm[0][1]->Fill(phi_pip, weight);
      phi_pim_cm[0][1]->Fill(phi_pim, weight);

      alpha_P_cm[0][1]->Fill(alpha_p, weight);
      alpha_pip_cm[0][1]->Fill(alpha_pip, weight);
      alpha_pim_cm[0][1]->Fill(alpha_pim, weight);
      // } else if (q2 > 5.0 && q2 < 15.0) {
      //   theta_P_cm[0][2]->Fill(theta_p, weight);
      //   theta_pip_cm[0][2]->Fill(theta_pip, weight);
      //   theta_pim_cm[0][2]->Fill(theta_pim, weight);
      //
      //   phi_P_cm[0][2]->Fill(phi_p, weight);
      //   phi_pip_cm[0][2]->Fill(phi_pip, weight);
      //   phi_pim_cm[0][2]->Fill(phi_pim, weight);
      //
      //   alpha_P_cm[0][2]->Fill(alpha_p, weight);
      //   alpha_pip_cm[0][2]->Fill(alpha_pip, weight);
      //   alpha_pim_cm[0][2]->Fill(alpha_pim, weight);
    }
  } else if (W > 2.5 && W < 3.5) {
    if (q2 < 4.5) {
      theta_P_cm[1][0]->Fill(theta_p, weight);
      theta_pip_cm[1][0]->Fill(theta_pip, weight);
      theta_pim_cm[1][0]->Fill(theta_pim, weight);

      phi_P_cm[1][0]->Fill(phi_p, weight);
      phi_pip_cm[1][0]->Fill(phi_pip, weight);
      phi_pim_cm[1][0]->Fill(phi_pim, weight);

      alpha_P_cm[1][0]->Fill(alpha_p, weight);
      alpha_pip_cm[1][0]->Fill(alpha_pip, weight);
      alpha_pim_cm[1][0]->Fill(alpha_pim, weight);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_cm[1][1]->Fill(theta_p, weight);
      theta_pip_cm[1][1]->Fill(theta_pip, weight);
      theta_pim_cm[1][1]->Fill(theta_pim, weight);

      phi_P_cm[1][1]->Fill(phi_p, weight);
      phi_pip_cm[1][1]->Fill(phi_pip, weight);
      phi_pim_cm[1][1]->Fill(phi_pim, weight);

      alpha_P_cm[1][1]->Fill(alpha_p, weight);
      alpha_pip_cm[1][1]->Fill(alpha_pip, weight);
      alpha_pim_cm[1][1]->Fill(alpha_pim, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     theta_P_cm[1][2]->Fill(theta_p, weight);
      //     theta_pip_cm[1][2]->Fill(theta_pip, weight);
      //     theta_pim_cm[1][2]->Fill(theta_pim, weight);
      //
      //     phi_P_cm[1][2]->Fill(phi_p, weight);
      //     phi_pip_cm[1][2]->Fill(phi_pip, weight);
      //     phi_pim_cm[1][2]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm[1][2]->Fill(alpha_p, weight);
      //     alpha_pip_cm[1][2]->Fill(alpha_pip, weight);
      //     alpha_pim_cm[1][2]->Fill(alpha_pim, weight);
      //   }
      // } else if (W > 2.5 && W < 3.0) {
      //   if (q2 < 2.0) {
      //     theta_P_cm[2][0]->Fill(theta_p, weight);
      //     theta_pip_cm[2][0]->Fill(theta_pip, weight);
      //     theta_pim_cm[2][0]->Fill(theta_pim, weight);
      //
      //     phi_P_cm[2][0]->Fill(phi_p, weight);
      //     phi_pip_cm[2][0]->Fill(phi_pip, weight);
      //     phi_pim_cm[2][0]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm[2][0]->Fill(alpha_p, weight);
      //     alpha_pip_cm[2][0]->Fill(alpha_pip, weight);
      //     alpha_pim_cm[2][0]->Fill(alpha_pim, weight);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_cm[2][1]->Fill(theta_p, weight);
      //     theta_pip_cm[2][1]->Fill(theta_pip, weight);
      //     theta_pim_cm[2][1]->Fill(theta_pim, weight);
      //
      //     phi_P_cm[2][1]->Fill(phi_p, weight);
      //     phi_pip_cm[2][1]->Fill(phi_pip, weight);
      //     phi_pim_cm[2][1]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm[2][1]->Fill(alpha_p, weight);
      //     alpha_pip_cm[2][1]->Fill(alpha_pip, weight);
      //     alpha_pim_cm[2][1]->Fill(alpha_pim, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     theta_P_cm[2][2]->Fill(theta_p, weight);
      //     theta_pip_cm[2][2]->Fill(theta_pip, weight);
      //     theta_pim_cm[2][2]->Fill(theta_pim, weight);
      //
      //     phi_P_cm[2][2]->Fill(phi_p, weight);
      //     phi_pip_cm[2][2]->Fill(phi_pip, weight);
      //     phi_pim_cm[2][2]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm[2][2]->Fill(alpha_p, weight);
      //     alpha_pip_cm[2][2]->Fill(alpha_pip, weight);
      //     alpha_pim_cm[2][2]->Fill(alpha_pim, weight);
      //   }
      // } else if (W > 3.0 && W < 3.5) {
      //   if (q2 < 2.0) {
      //     theta_P_cm[3][0]->Fill(theta_p, weight);
      //     theta_pip_cm[3][0]->Fill(theta_pip, weight);
      //     theta_pim_cm[3][0]->Fill(theta_pim, weight);
      //
      //     phi_P_cm[3][0]->Fill(phi_p, weight);
      //     phi_pip_cm[3][0]->Fill(phi_pip, weight);
      //     phi_pim_cm[3][0]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm[3][0]->Fill(alpha_p, weight);
      //     alpha_pip_cm[3][0]->Fill(alpha_pip, weight);
      //     alpha_pim_cm[3][0]->Fill(alpha_pim, weight);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_cm[3][1]->Fill(theta_p, weight);
      //     theta_pip_cm[3][1]->Fill(theta_pip, weight);
      //     theta_pim_cm[3][1]->Fill(theta_pim, weight);
      //
      //     phi_P_cm[3][1]->Fill(phi_p, weight);
      //     phi_pip_cm[3][1]->Fill(phi_pip, weight);
      //     phi_pim_cm[3][1]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm[3][1]->Fill(alpha_p, weight);
      //     alpha_pip_cm[3][1]->Fill(alpha_pip, weight);
      //     alpha_pim_cm[3][1]->Fill(alpha_pim, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     theta_P_cm[3][2]->Fill(theta_p, weight);
      //     theta_pip_cm[3][2]->Fill(theta_pip, weight);
      //     theta_pim_cm[3][2]->Fill(theta_pim, weight);
      //
      //     phi_P_cm[3][2]->Fill(phi_p, weight);
      //     phi_pip_cm[3][2]->Fill(phi_pip, weight);
      //     phi_pim_cm[3][2]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm[3][2]->Fill(alpha_p, weight);
      //     alpha_pip_cm[3][2]->Fill(alpha_pip, weight);
      //     alpha_pim_cm[3][2]->Fill(alpha_pim, weight);
    }
  }
}

void Histogram::Fill_Angles_thrown(float W, float q2, float theta_p, float theta_pip, float theta_pim, float phi_p,
                                   float phi_pip, float phi_pim, float alpha_p, float alpha_pip, float alpha_pim,
                                   float weight) {
  if (W <= 2.5) {
    if (q2 < 4.5) {
      theta_P_cm_thrown[0][0]->Fill(theta_p, weight);
      theta_pip_cm_thrown[0][0]->Fill(theta_pip, weight);
      theta_pim_cm_thrown[0][0]->Fill(theta_pim, weight);

      phi_P_cm_thrown[0][0]->Fill(phi_p, weight);
      phi_pip_cm_thrown[0][0]->Fill(phi_pip, weight);
      phi_pim_cm_thrown[0][0]->Fill(phi_pim, weight);

      alpha_P_cm_thrown[0][0]->Fill(alpha_p, weight);
      alpha_pip_cm_thrown[0][0]->Fill(alpha_pip, weight);
      alpha_pim_cm_thrown[0][0]->Fill(alpha_pim, weight);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_cm_thrown[0][1]->Fill(theta_p, weight);
      theta_pip_cm_thrown[0][1]->Fill(theta_pip, weight);
      theta_pim_cm_thrown[0][1]->Fill(theta_pim, weight);

      phi_P_cm_thrown[0][1]->Fill(phi_p, weight);
      phi_pip_cm_thrown[0][1]->Fill(phi_pip, weight);
      phi_pim_cm_thrown[0][1]->Fill(phi_pim, weight);

      alpha_P_cm_thrown[0][1]->Fill(alpha_p, weight);
      alpha_pip_cm_thrown[0][1]->Fill(alpha_pip, weight);
      alpha_pim_cm_thrown[0][1]->Fill(alpha_pim, weight);
      // } else if (q2 > 5.0 && q2 < 15.0) {
      //   theta_P_cm_thrown[0][2]->Fill(theta_p, weight);
      //   theta_pip_cm_thrown[0][2]->Fill(theta_pip, weight);
      //   theta_pim_cm_thrown[0][2]->Fill(theta_pim, weight);
      //
      //   phi_P_cm_thrown[0][2]->Fill(phi_p, weight);
      //   phi_pip_cm_thrown[0][2]->Fill(phi_pip, weight);
      //   phi_pim_cm_thrown[0][2]->Fill(phi_pim, weight);
      //
      //   alpha_P_cm_thrown[0][2]->Fill(alpha_p, weight);
      //   alpha_pip_cm_thrown[0][2]->Fill(alpha_pip, weight);
      //   alpha_pim_cm_thrown[0][2]->Fill(alpha_pim, weight);
    }
  } else if (W > 2.5 && W < 3.5) {
    if (q2 < 4.5) {
      theta_P_cm_thrown[1][0]->Fill(theta_p, weight);
      theta_pip_cm_thrown[1][0]->Fill(theta_pip, weight);
      theta_pim_cm_thrown[1][0]->Fill(theta_pim, weight);

      phi_P_cm_thrown[1][0]->Fill(phi_p, weight);
      phi_pip_cm_thrown[1][0]->Fill(phi_pip, weight);
      phi_pim_cm_thrown[1][0]->Fill(phi_pim, weight);

      alpha_P_cm_thrown[1][0]->Fill(alpha_p, weight);
      alpha_pip_cm_thrown[1][0]->Fill(alpha_pip, weight);
      alpha_pim_cm_thrown[1][0]->Fill(alpha_pim, weight);
    } else if (q2 > 4.5 && q2 < 15.0) {
      theta_P_cm_thrown[1][1]->Fill(theta_p, weight);
      theta_pip_cm_thrown[1][1]->Fill(theta_pip, weight);
      theta_pim_cm_thrown[1][1]->Fill(theta_pim, weight);

      phi_P_cm_thrown[1][1]->Fill(phi_p, weight);
      phi_pip_cm_thrown[1][1]->Fill(phi_pip, weight);
      phi_pim_cm_thrown[1][1]->Fill(phi_pim, weight);

      alpha_P_cm_thrown[1][1]->Fill(alpha_p, weight);
      alpha_pip_cm_thrown[1][1]->Fill(alpha_pip, weight);
      alpha_pim_cm_thrown[1][1]->Fill(alpha_pim, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     theta_P_cm_thrown[1][2]->Fill(theta_p, weight);
      //     theta_pip_cm_thrown[1][2]->Fill(theta_pip, weight);
      //     theta_pim_cm_thrown[1][2]->Fill(theta_pim, weight);
      //
      //     phi_P_cm_thrown[1][2]->Fill(phi_p, weight);
      //     phi_pip_cm_thrown[1][2]->Fill(phi_pip, weight);
      //     phi_pim_cm_thrown[1][2]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm_thrown[1][2]->Fill(alpha_p, weight);
      //     alpha_pip_cm_thrown[1][2]->Fill(alpha_pip, weight);
      //     alpha_pim_cm_thrown[1][2]->Fill(alpha_pim, weight);
      //   }
      // } else if (W > 2.5 && W < 3.0) {
      //   if (q2 < 2.0) {
      //     theta_P_cm_thrown[2][0]->Fill(theta_p, weight);
      //     theta_pip_cm_thrown[2][0]->Fill(theta_pip, weight);
      //     theta_pim_cm_thrown[2][0]->Fill(theta_pim, weight);
      //
      //     phi_P_cm_thrown[2][0]->Fill(phi_p, weight);
      //     phi_pip_cm_thrown[2][0]->Fill(phi_pip, weight);
      //     phi_pim_cm_thrown[2][0]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm_thrown[2][0]->Fill(alpha_p, weight);
      //     alpha_pip_cm_thrown[2][0]->Fill(alpha_pip, weight);
      //     alpha_pim_cm_thrown[2][0]->Fill(alpha_pim, weight);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_cm_thrown[2][1]->Fill(theta_p, weight);
      //     theta_pip_cm_thrown[2][1]->Fill(theta_pip, weight);
      //     theta_pim_cm_thrown[2][1]->Fill(theta_pim, weight);
      //
      //     phi_P_cm_thrown[2][1]->Fill(phi_p, weight);
      //     phi_pip_cm_thrown[2][1]->Fill(phi_pip, weight);
      //     phi_pim_cm_thrown[2][1]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm_thrown[2][1]->Fill(alpha_p, weight);
      //     alpha_pip_cm_thrown[2][1]->Fill(alpha_pip, weight);
      //     alpha_pim_cm_thrown[2][1]->Fill(alpha_pim, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     theta_P_cm_thrown[2][2]->Fill(theta_p, weight);
      //     theta_pip_cm_thrown[2][2]->Fill(theta_pip, weight);
      //     theta_pim_cm_thrown[2][2]->Fill(theta_pim, weight);
      //
      //     phi_P_cm_thrown[2][2]->Fill(phi_p, weight);
      //     phi_pip_cm_thrown[2][2]->Fill(phi_pip, weight);
      //     phi_pim_cm_thrown[2][2]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm_thrown[2][2]->Fill(alpha_p, weight);
      //     alpha_pip_cm_thrown[2][2]->Fill(alpha_pip, weight);
      //     alpha_pim_cm_thrown[2][2]->Fill(alpha_pim, weight);
      //   }
      // } else if (W > 3.0 && W < 3.5) {
      //   if (q2 < 2.0) {
      //     theta_P_cm_thrown[3][0]->Fill(theta_p, weight);
      //     theta_pip_cm_thrown[3][0]->Fill(theta_pip, weight);
      //     theta_pim_cm_thrown[3][0]->Fill(theta_pim, weight);
      //
      //     phi_P_cm_thrown[3][0]->Fill(phi_p, weight);
      //     phi_pip_cm_thrown[3][0]->Fill(phi_pip, weight);
      //     phi_pim_cm_thrown[3][0]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm_thrown[3][0]->Fill(alpha_p, weight);
      //     alpha_pip_cm_thrown[3][0]->Fill(alpha_pip, weight);
      //     alpha_pim_cm_thrown[3][0]->Fill(alpha_pim, weight);
      //   } else if (q2 > 2.0 && q2 < 5.0) {
      //     theta_P_cm_thrown[3][1]->Fill(theta_p, weight);
      //     theta_pip_cm_thrown[3][1]->Fill(theta_pip, weight);
      //     theta_pim_cm_thrown[3][1]->Fill(theta_pim, weight);
      //
      //     phi_P_cm_thrown[3][1]->Fill(phi_p, weight);
      //     phi_pip_cm_thrown[3][1]->Fill(phi_pip, weight);
      //     phi_pim_cm_thrown[3][1]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm_thrown[3][1]->Fill(alpha_p, weight);
      //     alpha_pip_cm_thrown[3][1]->Fill(alpha_pip, weight);
      //     alpha_pim_cm_thrown[3][1]->Fill(alpha_pim, weight);
      //   } else if (q2 > 5.0 && q2 < 15.0) {
      //     theta_P_cm_thrown[3][2]->Fill(theta_p, weight);
      //     theta_pip_cm_thrown[3][2]->Fill(theta_pip, weight);
      //     theta_pim_cm_thrown[3][2]->Fill(theta_pim, weight);
      //
      //     phi_P_cm_thrown[3][2]->Fill(phi_p, weight);
      //     phi_pip_cm_thrown[3][2]->Fill(phi_pip, weight);
      //     phi_pim_cm_thrown[3][2]->Fill(phi_pim, weight);
      //
      //     alpha_P_cm_thrown[3][2]->Fill(alpha_p, weight);
      //     alpha_pip_cm_thrown[3][2]->Fill(alpha_pip, weight);
      //     alpha_pim_cm_thrown[3][2]->Fill(alpha_pim, weight);
    }
  }
}

void Histogram::Fill_W_hist_Xpip_all_sec(float W, float weight) { W_hist_Xpip_all_sec->Fill(W, weight); }
void Histogram::Fill_WvsmmSQ_2pi(float W, float W_dpp, float delta_zero_, float rho_, float mmSQ, int sec_number,
                                 float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_hist_2pi[sec_number]->Fill(W, wt);
    W_hist_delta_pp[sec_number]->Fill(W_dpp, wt);
    W_hist_delta_zero[sec_number]->Fill(delta_zero_, wt);
    W_hist_rho[sec_number]->Fill(rho_, wt);
    W_vs_mmSQ_2pi[sec_number][0]->Fill(W, mmSQ, wt);
  }
}
void Histogram::Fill_WvsmmSQ_singlepip(float W, float mmSQ, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_hist_singlepip[sec_number]->Fill(W, wt);
    W_vs_mmSQ_singlepip[sec_number][0]->Fill(W, mmSQ, wt);
  }
}
void Histogram::Fill_WvsmmSQ_anti_ep(float W, float mmSQ, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_vs_mmSQ_ep[sec_number][1]->Fill(W, mmSQ, wt);
  }
}
void Histogram::Fill_WvsmmSQ_anti_2pi(float W, float W_dpp, float delta_zero_, float rho_, float mmSQ, int sec_number,
                                      float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_vs_mmSQ_2pi[sec_number][1]->Fill(W, mmSQ, wt);
  }
}
void Histogram::Fill_WvsmmSQ_anti_singlepip(float W, float mmSQ, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_vs_mmSQ_singlepip[sec_number][1]->Fill(W, mmSQ, wt);
  }
}

// void Histogram::Make_hist_cc() {
//   for (int i = 0; i < cc_num; i++) {
//     hname.append("cc_total_");
//     htitle.append("cc_total_");
//     hname.append(cc_name[i]);
//     htitle.append(cc_name[i]);
//     cherenkov_total[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0,
//     60); hname.clear(); htitle.clear();
//
//     hname.append("cc_ltcc_");
//     htitle.append("cc_ltcc_");
//     hname.append(cc_name[i]);
//     htitle.append(cc_name[i]);
//     cherenkov_ltcc[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 20);
//     hname.clear();
//     htitle.clear();
//
//     hname.append("cc_htcc_");
//     htitle.append("cc_htcc_");
//     hname.append(cc_name[i]);
//     htitle.append(cc_name[i]);
//     cherenkov_htcc[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 60);
//     hname.clear();
//     htitle.clear();
//   }
// }
// void Histogram::Fill_hist_cc_tot(float tot_el) {
// cherenkov_total[0]->Fill(tot_el); } void Histogram::Fill_hist_cc_ltcc(float
// ltcc_el) { cherenkov_ltcc[0]->Fill(ltcc_el); } void
// Histogram::Fill_hist_cc_htcc(float htcc_el) {
// cherenkov_htcc[0]->Fill(htcc_el); } void
// Histogram::Fill_hist_cc_tot_pim(float tot_pim) {
// cherenkov_total[1]->Fill(tot_pim); } void
// Histogram::Fill_hist_cc_ltcc_pim(float ltcc_pim) {
// cherenkov_ltcc[1]->Fill(ltcc_pim); } void
// Histogram::Fill_hist_cc_htcc_pim(float htcc_pim) {
// cherenkov_htcc[1]->Fill(htcc_pim); } void
// Histogram::Fill_hist_cc_tot_pip(float tot_pip) {
// cherenkov_total[2]->Fill(tot_pip); } void
// Histogram::Fill_hist_cc_ltcc_pip(float ltcc_pip) {
// cherenkov_ltcc[2]->Fill(ltcc_pip); } void
// Histogram::Fill_hist_cc_htcc_pip(float htcc_pip) {
// cherenkov_htcc[2]->Fill(htcc_pip); }
//  } else if (cc_num == 1) {

// void Histogram::Write_hist_cc() {
//   for (int i = 0; i < cc_num; i++) {
//     cherenkov_total[i]->SetXTitle("cc total");
//     cherenkov_total[i]->Write();
//     delete cherenkov_total[i];
//
//     cherenkov_ltcc[i]->SetXTitle("cc ltcc");
//     cherenkov_ltcc[i]->Write();
//     delete cherenkov_ltcc[i];
//
//     cherenkov_htcc[i]->SetXTitle("cc htcc");
//     cherenkov_htcc[i]->Write();
//     delete cherenkov_htcc[i];
//   }
// }

void Histogram::Fill_lu_dist(float li) { lu_side_distribution->Fill(li); }
void Histogram::Fill_lv_dist(float li) { lv_side_distribution->Fill(li); }
void Histogram::Fill_lw_dist(float li) { lw_side_distribution->Fill(li); }

void Histogram::Fill_vertex_vz(float vz) { vertex_vz->Fill(vz); }

void Histogram::Fill_pid_size_fn(int pid_size, int pid) {
  pid_size1->Fill(pid_size);
  pid_check->Fill(pid);
}
void Histogram::Fill_mom_diff_e_pim(float mom_diff, float mom_diff_x, float mom_diff_y, float mom_diff_z) {
  mom_diff_e_pim->Fill(mom_diff);
  mom_diff_e_pim_x->Fill(mom_diff_x);
  mom_diff_e_pim_y->Fill(mom_diff_y);
  mom_diff_e_pim_z->Fill(mom_diff_z);
}
void Histogram::Fill_mom_diff_p_pip(float mom_diff, float mom_diff_x, float mom_diff_y, float mom_diff_z) {
  mom_diff_p_pip->Fill(mom_diff);
  mom_diff_p_pip_x->Fill(mom_diff_x);
  mom_diff_p_pip_y->Fill(mom_diff_y);
  mom_diff_p_pip_z->Fill(mom_diff_z);
}
void Histogram::Fill_mom_ele(float mom, float mom_z) {
  mom_e->Fill(mom);
  mom_e_z->Fill(mom_z);
}
void Histogram::Fill_mom_p(float mom, float mom_z) {
  mom_p->Fill(mom);
  mom_p_z->Fill(mom_z);
}
void Histogram::Fill_mom_pip(float mom, float mom_z) {
  mom_pip->Fill(mom);
  mom_pip_z->Fill(mom_z);
}
void Histogram::Fill_mom_pim(float mom, float mom_z) {
  mom_pim->Fill(mom);
  mom_pim_z->Fill(mom_z);
}
// void Histogram::Fill_theta_P(float theta_p, float theta_pip_, float
// theta_pim_) {
//         theta_prot->Fill(theta_p);
//         theta_pip->Fill(theta_pip_);
//         theta_pim->Fill(theta_pim_);
// }
// void Histogram::Fill_Phi_cm(float Phi_p, float Phi_pip_, float Phi_pim_) {
//         Phi_prot->Fill(Phi_p);
//         Phi_pip->Fill(Phi_pip_);
//         Phi_pim->Fill(Phi_pim_);
// }

void Histogram::Fill_W_vs_Q2_thrown(float w, float q2, float wt) {
  W_thrown->Fill(w, wt);
  W_vs_Q2_thrown->Fill(w, q2, wt);
}
void Histogram::Fill_W_vs_Q2_all_sec(float W, float q2, float wt) {
  // if (W <= 2.0)
  W_vs_Q2_all_sec->Fill(W, q2, wt);

  // else if (W >= 2.0 && W <= 2.5)
  //   W_vs_Q2_all_sec[1]->Fill(W, q2, wt);
  //
  // else if (W >= 2.5 && W <= 3.0)
  //   W_vs_Q2_all_sec[2]->Fill(W, q2, wt);
  //
  // else if (W >= 3.0 && W <= 3.5)
  //   W_vs_Q2_all_sec[3]->Fill(W, q2, wt);
}
// void Histogram::Fill_hist_mass_vs_q2_prot(float W, float m_p, float theta, float q2, float wt) {
//   if (W <= 2.0) {
//     if (q2 < 2.0) {
//       Prot_theta_lab_vs_W[0][0]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[0][0]->Fill(m_p, q2, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Prot_theta_lab_vs_W[0][1]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[0][1]->Fill(m_p, q2, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Prot_theta_lab_vs_W[0][2]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[0][2]->Fill(m_p, q2, wt);
//     }
//   } else if (W > 2.0 && W < 2.5) {
//     if (q2 < 2.0) {
//       Prot_theta_lab_vs_W[1][0]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[1][0]->Fill(m_p, q2, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Prot_theta_lab_vs_W[1][1]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[1][1]->Fill(m_p, q2, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Prot_theta_lab_vs_W[1][2]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[1][2]->Fill(m_p, q2, wt);
//     }
//
//   } else if (W >= 2.5 && W <= 3.0) {
//     if (q2 < 2.0) {
//       Prot_theta_lab_vs_W[2][0]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[2][0]->Fill(m_p, q2, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Prot_theta_lab_vs_W[2][1]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[2][1]->Fill(m_p, q2, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Prot_theta_lab_vs_W[2][2]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[2][2]->Fill(m_p, q2, wt);
//     }
//   } else if (W >= 3.0 && W <= 3.5) {
//     if (q2 < 2.0) {
//       Prot_theta_lab_vs_W[3][0]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[3][0]->Fill(m_p, q2, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Prot_theta_lab_vs_W[3][1]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[3][1]->Fill(m_p, q2, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Prot_theta_lab_vs_W[3][2]->Fill(W, theta, wt);
//       Prot_mass_w_vs_Q2[3][2]->Fill(m_p, q2, wt);
//     }
//   }
// }
// void Histogram::Fill_hist_mass_vs_q2_pip(float W, float m_pip, float theta, float q2, float wt) {
//   if (W <= 2.0) {
//     if (q2 < 2.0) {
//       Pip_theta_lab_vs_W[0][0]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[0][0]->Fill(m_pip, q2, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Pip_theta_lab_vs_W[0][1]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[0][1]->Fill(m_pip, q2, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Pip_theta_lab_vs_W[0][2]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[0][2]->Fill(m_pip, q2, wt);
//     }
//   } else if (W > 2.0 && W < 2.5) {
//     if (q2 < 2.0) {
//       Pip_theta_lab_vs_W[1][0]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[1][0]->Fill(m_pip, q2, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Pip_theta_lab_vs_W[1][1]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[1][1]->Fill(m_pip, q2, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Pip_theta_lab_vs_W[1][2]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[1][2]->Fill(m_pip, q2, wt);
//     }
//
//   } else if (W >= 2.5 && W <= 3.0) {
//     if (q2 < 2.0) {
//       Pip_theta_lab_vs_W[2][0]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[2][0]->Fill(m_pip, q2, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Pip_theta_lab_vs_W[2][1]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[2][1]->Fill(m_pip, q2, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Pip_theta_lab_vs_W[2][2]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[2][2]->Fill(m_pip, q2, wt);
//     }
//   } else if (W >= 3.0 && W <= 3.5) {
//     if (q2 < 2.0) {
//       Pip_theta_lab_vs_W[3][0]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[3][0]->Fill(m_pip, q2, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Pip_theta_lab_vs_W[3][1]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[3][1]->Fill(m_pip, q2, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Pip_theta_lab_vs_W[3][2]->Fill(W, theta, wt);
//       Pip_mass_w_vs_Q2[3][2]->Fill(m_pip, q2, wt);
//     }
//   }
// }
// void Histogram::Fill_hist_mass_vs_q2_pim(float W, float m_pim, float theta, float q2, float wt) {
//   if (W <= 2.0) {
//     if (q2 < 2.0) {
//       Pim_theta_lab_vs_W[0][0]->Fill(W, theta, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Pim_theta_lab_vs_W[0][1]->Fill(W, theta, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Pim_theta_lab_vs_W[0][2]->Fill(W, theta, wt);
//     }
//   } else if (W > 2.0 && W < 2.5) {
//     if (q2 < 2.0) {
//       Pim_theta_lab_vs_W[1][0]->Fill(W, theta, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Pim_theta_lab_vs_W[1][1]->Fill(W, theta, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Pim_theta_lab_vs_W[1][2]->Fill(W, theta, wt);
//     }
//
//   } else if (W >= 2.5 && W <= 3.0) {
//     if (q2 < 2.0) {
//       Pim_theta_lab_vs_W[2][0]->Fill(W, theta, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Pim_theta_lab_vs_W[2][1]->Fill(W, theta, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Pim_theta_lab_vs_W[2][2]->Fill(W, theta, wt);
//     }
//   } else if (W >= 3.0 && W <= 3.5) {
//     if (q2 < 2.0) {
//       Pim_theta_lab_vs_W[3][0]->Fill(W, theta, wt);
//     } else if (q2 > 2.0 && q2 < 5.0) {
//       Pim_theta_lab_vs_W[3][1]->Fill(W, theta, wt);
//     } else if (q2 > 5.0 && q2 < 10.0) {
//       Pim_theta_lab_vs_W[3][2]->Fill(W, theta, wt);
//     }
//   }
// }
void Histogram::Write_Angles() {
  for (int i = 0; i < w_range_num; i++) {
    for (int j = 0; j < q2_range_num; j++) {
      theta_P_cm_thrown[i][j]->SetXTitle("theta_P_cm_thrown (Deg)");
      theta_P_cm_thrown[i][j]->Write();
      delete theta_P_cm_thrown[i][j];
      phi_P_cm_thrown[i][j]->SetXTitle("phi_P_cm_thrown (Deg)");
      phi_P_cm_thrown[i][j]->Write();
      delete phi_P_cm_thrown[i][j];
      alpha_P_cm_thrown[i][j]->SetXTitle("alpha_P_cm_thrown (Deg)");
      alpha_P_cm_thrown[i][j]->Write();
      delete alpha_P_cm_thrown[i][j];

      theta_P_cm[i][j]->SetXTitle("theta_P_cm (Deg)");
      theta_P_cm[i][j]->Write();
      delete theta_P_cm[i][j];
      phi_P_cm[i][j]->SetXTitle("phi_P_cm (Deg)");
      phi_P_cm[i][j]->Write();
      delete phi_P_cm[i][j];
      alpha_P_cm[i][j]->SetXTitle("alpha_P_cm (Deg)");
      alpha_P_cm[i][j]->Write();
      delete alpha_P_cm[i][j];

      theta_pip_cm_thrown[i][j]->SetXTitle("theta_pip_cm_thrown (Deg)");
      theta_pip_cm_thrown[i][j]->Write();
      delete theta_pip_cm_thrown[i][j];
      phi_pip_cm_thrown[i][j]->SetXTitle("phi_pip_cm_thrown (Deg)");
      phi_pip_cm_thrown[i][j]->Write();
      delete phi_pip_cm_thrown[i][j];
      alpha_pip_cm_thrown[i][j]->SetXTitle("alpha_pip_cm_thrown (Deg)");
      alpha_pip_cm_thrown[i][j]->Write();
      delete alpha_pip_cm_thrown[i][j];

      theta_pip_cm[i][j]->SetXTitle("theta_pip_cm (Deg)");
      theta_pip_cm[i][j]->Write();
      delete theta_pip_cm[i][j];
      phi_pip_cm[i][j]->SetXTitle("phi_pip_cm (Deg)");
      phi_pip_cm[i][j]->Write();
      delete phi_pip_cm[i][j];
      alpha_pip_cm[i][j]->SetXTitle("alpha_pip_cm (Deg)");
      alpha_pip_cm[i][j]->Write();
      delete alpha_pip_cm[i][j];

      theta_pim_cm_thrown[i][j]->SetXTitle("theta_pim_cm_thrown (Deg)");
      theta_pim_cm_thrown[i][j]->Write();
      delete theta_pim_cm_thrown[i][j];
      phi_pim_cm_thrown[i][j]->SetXTitle("phi_pim_cm_thrown (Deg)");
      phi_pim_cm_thrown[i][j]->Write();
      delete phi_pim_cm_thrown[i][j];
      alpha_pim_cm_thrown[i][j]->SetXTitle("alpha_pim_cm_thrown (Deg)");
      alpha_pim_cm_thrown[i][j]->Write();
      delete alpha_pim_cm_thrown[i][j];

      theta_pim_cm[i][j]->SetXTitle("theta_pim_cm (Deg)");
      theta_pim_cm[i][j]->Write();
      delete theta_pim_cm[i][j];
      phi_pim_cm[i][j]->SetXTitle("phi_pim_cm (Deg)");
      phi_pim_cm[i][j]->Write();
      delete phi_pim_cm[i][j];
      alpha_pim_cm[i][j]->SetXTitle("alpha_pim_cm (Deg)");
      alpha_pim_cm[i][j]->Write();
      delete alpha_pim_cm[i][j];
    }
  }
}
void Histogram::Write_WvsQ2() {
  W_thrown->SetXTitle("W_thrown (GeV)");
  W_thrown->Write();
  W_vs_Q2_thrown->SetXTitle("W_thrown (GeV)");
  W_vs_Q2_thrown->SetYTitle("Q2");
  W_vs_Q2_thrown->SetOption("COLZ");
  W_vs_Q2_thrown->Write();
  W_vs_Q2_all_sec->SetXTitle("W (GeV)");
  W_vs_Q2_all_sec->SetYTitle("Q2");
  W_vs_Q2_all_sec->SetOption("COLZ");
  W_vs_Q2_all_sec->Write();
  // delete W_vs_Q2_all_sec[i];
  for (int i = 0; i < w_range_num; i++) {
    for (int j = 0; j < q2_range_num; j++) {
      Prot_mass_w_vs_Q2[i][j]->SetXTitle("prot_mass (GeV)");
      Prot_mass_w_vs_Q2[i][j]->SetYTitle("Q2");
      Prot_mass_w_vs_Q2[i][j]->SetOption("COLZ");
      Prot_mass_w_vs_Q2[i][j]->Write();
      delete Prot_mass_w_vs_Q2[i][j];
      Pip_mass_w_vs_Q2[i][j]->SetXTitle("pip_mass (GeV)");
      Pip_mass_w_vs_Q2[i][j]->SetYTitle("Q2");
      Pip_mass_w_vs_Q2[i][j]->SetOption("COLZ");
      Pip_mass_w_vs_Q2[i][j]->Write();
      delete Pip_mass_w_vs_Q2[i][j];
      Prot_theta_lab_vs_W[i][j]->SetXTitle("W (GeV)");
      Prot_theta_lab_vs_W[i][j]->SetYTitle("prot_theta_lab");
      Prot_theta_lab_vs_W[i][j]->SetOption("COLZ");
      Prot_theta_lab_vs_W[i][j]->Write();
      delete Prot_theta_lab_vs_W[i][j];
      Pip_theta_lab_vs_W[i][j]->SetXTitle("W (GeV)");
      Pip_theta_lab_vs_W[i][j]->SetYTitle("pip_theta_lab");
      Pip_theta_lab_vs_W[i][j]->SetOption("COLZ");
      Pip_theta_lab_vs_W[i][j]->Write();
      delete Pip_theta_lab_vs_W[i][j];
      Pim_theta_lab_vs_W[i][j]->SetXTitle("W (GeV)");
      Pim_theta_lab_vs_W[i][j]->SetYTitle("pim_theta_lab");
      Pim_theta_lab_vs_W[i][j]->SetOption("COLZ");
      Pim_theta_lab_vs_W[i][j]->Write();
      delete Pim_theta_lab_vs_W[i][j];

      W_hist_twopi_all_sec[i][j]->SetXTitle("W_P_pip_pim (GeV)");
      W_hist_twopi_all_sec[i][j]->Write();
      delete W_hist_twopi_all_sec[i][j];
      inv_mass_P_pip_all_sec[i][j]->SetXTitle("inv_mass_P_pip (GeV)");
      inv_mass_P_pip_all_sec[i][j]->Write();
      delete inv_mass_P_pip_all_sec[i][j];
      inv_mass_P_pim_all_sec[i][j]->SetXTitle("inv_mass_P_pim (GeV)");
      inv_mass_P_pim_all_sec[i][j]->Write();
      delete inv_mass_P_pim_all_sec[i][j];
      inv_mass_pip_pim_all_sec[i][j]->SetXTitle("inv_mass_pip_pim (GeV)");
      inv_mass_pip_pim_all_sec[i][j]->Write();
      delete inv_mass_pip_pim_all_sec[i][j];

      theta_pim_vs_mass_Ppip[i][j]->SetXTitle("Inv_m_P_pi+");
      theta_pim_vs_mass_Ppip[i][j]->SetYTitle("theta_pi_cm");
      theta_pim_vs_mass_Ppip[i][j]->SetOption("COLZ");
      theta_pim_vs_mass_Ppip[i][j]->Write();
      delete theta_pim_vs_mass_Ppip[i][j];

      theta_pip_vs_mass_Ppim[i][j]->SetXTitle("Inv_m_P-pi-");
      theta_pip_vs_mass_Ppim[i][j]->SetYTitle("theta_pi+_cm");
      theta_pip_vs_mass_Ppim[i][j]->SetOption("COLZ");
      theta_pip_vs_mass_Ppim[i][j]->Write();
      delete theta_pip_vs_mass_Ppim[i][j];

      theta_P_vs_mass_pip_pim[i][j]->SetXTitle("Inv_m_Pi+pi-");
      theta_P_vs_mass_pip_pim[i][j]->SetYTitle("theta_P_cm");
      theta_P_vs_mass_pip_pim[i][j]->SetOption("COLZ");
      theta_P_vs_mass_pip_pim[i][j]->Write();
      delete theta_P_vs_mass_pip_pim[i][j];

      theta_pim_lab_vs_mass_Ppip[i][j]->SetXTitle("Inv_m_P_pi+");
      theta_pim_lab_vs_mass_Ppip[i][j]->SetYTitle("theta_pi-_lab");
      theta_pim_lab_vs_mass_Ppip[i][j]->SetOption("COLZ");
      theta_pim_lab_vs_mass_Ppip[i][j]->Write();
      delete theta_pim_lab_vs_mass_Ppip[i][j];

      theta_pip_lab_vs_mass_Ppim[i][j]->SetXTitle("Inv_m_P-pi-");
      theta_pip_lab_vs_mass_Ppim[i][j]->SetYTitle("theta_pi+_lab");
      theta_pip_lab_vs_mass_Ppim[i][j]->SetOption("COLZ");
      theta_pip_lab_vs_mass_Ppim[i][j]->Write();
      delete theta_pip_lab_vs_mass_Ppim[i][j];

      theta_P_lab_vs_mass_pip_pim[i][j]->SetXTitle("Inv_m_Pi+pi-");
      theta_P_lab_vs_mass_pip_pim[i][j]->SetYTitle("theta_P_lab");
      theta_P_lab_vs_mass_pip_pim[i][j]->SetOption("COLZ");
      theta_P_lab_vs_mass_pip_pim[i][j]->Write();
      delete theta_P_lab_vs_mass_pip_pim[i][j];

      W_hist_twopi_thrown[i][j]->SetXTitle("W_P_pip_pim_thrown (GeV)");
      W_hist_twopi_thrown[i][j]->Write();
      delete W_hist_twopi_thrown[i][j];
      inv_mass_P_pip_thrown[i][j]->SetXTitle("inv_mass_P_pip_thrown (GeV)");
      inv_mass_P_pip_thrown[i][j]->Write();
      delete inv_mass_P_pip_thrown[i][j];
      inv_mass_P_pim_thrown[i][j]->SetXTitle("inv_mass_P_pim_thrown (GeV)");
      inv_mass_P_pim_thrown[i][j]->Write();
      delete inv_mass_P_pim_thrown[i][j];
      inv_mass_pip_pim_thrown[i][j]->SetXTitle("inv_mass_pip_pim_thrown (GeV)");
      inv_mass_pip_pim_thrown[i][j]->Write();
      delete inv_mass_pip_pim_thrown[i][j];

      theta_pim_vs_mass_Ppip_thrown[i][j]->SetXTitle("Inv_m_P_pi+");
      theta_pim_vs_mass_Ppip_thrown[i][j]->SetYTitle("theta_pi_cm");
      theta_pim_vs_mass_Ppip_thrown[i][j]->SetOption("COLZ");
      theta_pim_vs_mass_Ppip_thrown[i][j]->Write();
      delete theta_pim_vs_mass_Ppip_thrown[i][j];

      theta_pip_vs_mass_Ppim_thrown[i][j]->SetXTitle("Inv_m_P-pi-");
      theta_pip_vs_mass_Ppim_thrown[i][j]->SetYTitle("theta_pi+_cm");
      theta_pip_vs_mass_Ppim_thrown[i][j]->SetOption("COLZ");
      theta_pip_vs_mass_Ppim_thrown[i][j]->Write();
      delete theta_pip_vs_mass_Ppim_thrown[i][j];

      theta_P_vs_mass_pip_pim_thrown[i][j]->SetXTitle("Inv_m_Pi+pi-");
      theta_P_vs_mass_pip_pim_thrown[i][j]->SetYTitle("theta_P_cm");
      theta_P_vs_mass_pip_pim_thrown[i][j]->SetOption("COLZ");
      theta_P_vs_mass_pip_pim_thrown[i][j]->Write();
      delete theta_P_vs_mass_pip_pim_thrown[i][j];

      theta_pim_lab_vs_mass_Ppip_thrown[i][j]->SetXTitle("Inv_m_P_pi+");
      theta_pim_lab_vs_mass_Ppip_thrown[i][j]->SetYTitle("theta_pi-_lab");
      theta_pim_lab_vs_mass_Ppip_thrown[i][j]->SetOption("COLZ");
      theta_pim_lab_vs_mass_Ppip_thrown[i][j]->Write();
      delete theta_pim_lab_vs_mass_Ppip_thrown[i][j];

      theta_pip_lab_vs_mass_Ppim_thrown[i][j]->SetXTitle("Inv_m_P-pi-");
      theta_pip_lab_vs_mass_Ppim_thrown[i][j]->SetYTitle("theta_pi+_lab");
      theta_pip_lab_vs_mass_Ppim_thrown[i][j]->SetOption("COLZ");
      theta_pip_lab_vs_mass_Ppim_thrown[i][j]->Write();
      delete theta_pip_lab_vs_mass_Ppim_thrown[i][j];

      theta_P_lab_vs_mass_pip_pim_thrown[i][j]->SetXTitle("Inv_m_Pi+pi-");
      theta_P_lab_vs_mass_pip_pim_thrown[i][j]->SetYTitle("theta_P_lab");
      theta_P_lab_vs_mass_pip_pim_thrown[i][j]->SetOption("COLZ");
      theta_P_lab_vs_mass_pip_pim_thrown[i][j]->Write();
      delete theta_P_lab_vs_mass_pip_pim_thrown[i][j];
    }
  }
  W_hist_Xpip_all_sec->SetXTitle("W_Xpip_all_sec (GeV)");
  W_hist_Xpip_all_sec->Write();
  // delete W_hist_Xpip_all_sec;

  for (int i = 0; i < sec_num; i++) {
    W_vs_Q2[i]->SetXTitle("W (GeV)");
    W_vs_Q2[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_Q2[i]->SetOption("COLZ");
    W_vs_Q2[i]->Write();
    delete W_vs_Q2[i];

    W_hist[i]->SetXTitle("W (GeV)");
    W_hist[i]->Write();
    delete W_hist[i];

    W_hist_ep[i]->SetXTitle("W_ep (GeV)");
    W_hist_ep[i]->Write();
    delete W_hist_ep[i];

    W_hist_2pi[i]->SetXTitle("W_2pi (GeV)");
    W_hist_2pi[i]->Write();
    delete W_hist_2pi[i];

    W_hist_delta_pp[i]->SetXTitle("inv_mass_p#pi+ (GeV)");
    W_hist_delta_pp[i]->Write();
    delete W_hist_delta_pp[i];

    W_hist_delta_zero[i]->SetXTitle("inv_mass_p#pi- (GeV)");
    W_hist_delta_zero[i]->Write();
    delete W_hist_delta_zero[i];

    W_hist_rho[i]->SetXTitle("inv_mass_#pi+p#pi- (GeV)");
    W_hist_rho[i]->Write();
    delete W_hist_rho[i];

    W_hist_singlepip[i]->SetXTitle("W_n#pi+ (GeV)");
    W_hist_singlepip[i]->Write();
    delete W_hist_singlepip[i];

    Q2_hist[i]->SetXTitle("Q^{2} (GeV^{2})");
    Q2_hist[i]->Write();
    delete Q2_hist[i];

    for (int j = 0; j < cut_y_n; j++) {
      W_vs_mmSQ_ep[i][j]->SetXTitle("W (GeV)");
      W_vs_mmSQ_ep[i][j]->SetYTitle("mm^{2} (GeV^{2})");
      W_vs_mmSQ_ep[i][j]->SetOption("COLZ");
      W_vs_mmSQ_ep[i][j]->Write();
      delete W_vs_mmSQ_ep[i][j];

      W_vs_mmSQ_2pi[i][j]->SetXTitle("W (GeV)");
      W_vs_mmSQ_2pi[i][j]->SetYTitle("mm^{2} (GeV^{2})");
      W_vs_mmSQ_2pi[i][j]->SetOption("COLZ");
      W_vs_mmSQ_2pi[i][j]->Write();
      delete W_vs_mmSQ_2pi[i][j];

      W_vs_mmSQ_singlepip[i][j]->SetXTitle("W (GeV)");
      W_vs_mmSQ_singlepip[i][j]->SetYTitle("mm^{2} (GeV^{2})");
      W_vs_mmSQ_singlepip[i][j]->SetOption("COLZ");
      W_vs_mmSQ_singlepip[i][j]->Write();
      delete W_vs_mmSQ_singlepip[i][j];
    }
  }
}
void Histogram::makeHists_MM() {
  for (size_t m = 0; m < mm_num; m++) {
    for (size_t e = 0; e < mm_events_num; e++) {
      for (int i = 0; i < sec_num; i++) {
        hname.append(mm_name[m]);
        htitle.append(mm_name[m]);
        hname.append("_");
        htitle.append(" ");
        hname.append(mm_events_name[e]);
        htitle.append(mm_events_name[e]);
        hname.append("_events_");
        htitle.append(" events ");
        hname.append(sec_name[i]);
        htitle.append(sec_name[i]);
        MM_hist[m][e][i] =
            new TH1D(hname.c_str(), htitle.c_str(), bins, Histogram::mm_lim_min(m, e), Histogram::mm_lim_max(m, e));
        hname.clear();
        htitle.clear();
      }
    }
  }
}

void Histogram::Fill_ep_mm(float mm, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][0][sec_number]->Fill(mm, wt);
  }
}
void Histogram::Fill_ep_mmSQ(float mm, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][0][sec_number]->Fill(mm, wt);
  }
}
void Histogram::Fill_2pion_mm(float mm, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][1][sec_number]->Fill(mm, wt);
  }
}
void Histogram::Fill_2pion_mmSQ(float mm, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][1][sec_number]->Fill(mm, wt);
  }
}
void Histogram::Fill_pip_mm(float mm, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][2][sec_number]->Fill(mm, wt);
  }
}
void Histogram::Fill_pip_mmSQ(float mm, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][2][sec_number]->Fill(mm, wt);
  }
}
void Histogram::Fill_pim_mm(float mm, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][3][sec_number]->Fill(mm, wt);
  }
}
void Histogram::Fill_pim_mmSQ(float mm, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][3][sec_number]->Fill(mm, wt);
  }
}
void Histogram::Fill_MM_wop_e_prime(float mm_1, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][4][sec_number]->Fill(mm_1, wt);
  }
}
void Histogram::Fill_MMSQ_wop_e_prime(float mm_1, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][4][sec_number]->Fill(mm_1, wt);
  }
}
void Histogram::Fill_MM_wop_2pion(float mm_1, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][5][sec_number]->Fill(mm_1, wt);
  }
}
void Histogram::Fill_MMSQ_wop_2pion(float mm_1, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][5][sec_number]->Fill(mm_1, wt);
  }
}
void Histogram::Fill_MM_wop_pip(float mm_1, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][6][sec_number]->Fill(mm_1, wt);
  }
}
void Histogram::Fill_MMSQ_wop_pip(float mm_1, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][6][sec_number]->Fill(mm_1, wt);
  }
}
void Histogram::Fill_MM_wop_pim(float mm_1, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][7][sec_number]->Fill(mm_1, wt);
  }
}
void Histogram::Fill_MMSQ_wop_pim(float mm_1, int sec_number, float wt) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][7][sec_number]->Fill(mm_1, wt);
  }
}

void Histogram::Write_MM_hist() {
  for (size_t m = 0; m < mm_num; m++) {
    for (size_t e = 0; e < mm_events_num; e++) {
      for (int i = 0; i < sec_num; i++) {
        // if (m == 1 && e == 1) {
        //   MM_hist[m][e][i]->Fit("gaus", "", "", -0.05, 0.05);
        // }
        MM_hist[m][e][i]->Write();
        delete MM_hist[m][e][i];
      }
    }
  }
}
// void Histogram::makeHists_deltat() {
//         for (size_t i = 0; i < with_id_num; i++) {
//                 hname.append("delta_t_vertex");
//                 htitle.append("#Deltat vertex particle");
//                 hname.append("_");
//                 htitle.append(" ");
//                 hname.append(id_name[i]);
//                 htitle.append(id_name[i]);
//                 delta_t_vertex[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, Dt_min,
//                 Dt_max); hname.clear(); htitle.clear();
//         }
//
//         for (size_t p = 0; p < particle_num; p++) {
//                 for (size_t c = 0; c < charge_num; c++) {
//                         for (size_t i = 0; i < with_id_num; i++) {
//                                 hname.append("delta_t_");
//                                 htitle.append("#Deltat ");
//                                 hname.append(particle_name[p]);
//                                 htitle.append(particle_name[p]);
//                                 hname.append("_");
//                                 htitle.append(" ");
//                                 hname.append(charge_name[c]);
//                                 htitle.append(charge_name[c]);
//                                 hname.append("_");
//                                 htitle.append(" ");
//                                 hname.append(id_name[i]);
//                                 htitle.append(id_name[i]);
//                                 delta_t_hist[p][c][i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min,
//                                 p_max, bins, Dt_min, Dt_max); hname.clear(); htitle.clear();
//                         }
//                 }
//         }
//         for (size_t p = 0; p < particle_num; p++) {
//                 for (size_t c = 0; c < charge_num; c++) {
//                         for (size_t ct = 0; ct < cut_without_cut_num; ct++) {
//                                 hname.append("delta_t_");
//                                 htitle.append("#Deltat ");
//                                 hname.append(particle_name[p]);
//                                 htitle.append(particle_name[p]);
//                                 hname.append("_ctof_");
//                                 htitle.append(" ctof ");
//                                 hname.append(charge_name[c]);
//                                 htitle.append(charge_name[c]);
//                                 hname.append("_");
//                                 htitle.append(" ");
//                                 hname.append(cut_without_cut_name[ct]);
//                                 htitle.append(cut_without_cut_name[ct]);
//                                 delta_t_hist_ctof[p][c][ct] = new TH2D(hname.c_str(), htitle.c_str(), bins,
//                                 p_min, p_max, bins, Dt_min, Dt_max); hname.clear(); htitle.clear();
//                         }
//                 }
//         }
// }
// void Histogram::Fill_deltat_prot_after(int pid, int charge, float dt_p, float dt_pip, float momentum) {
//         Deltat_P_after_mmsq_cut->Fill(momentum, dt_p);
//         Deltat_pip_after_mmsq_cut->Fill(momentum, dt_pip);
// }
// void Histogram::Fill_deltat_ctof_prot_after(int pid, int charge, float dt_ctof_p, float dt_ctof_pip, float
// momentum)
// {
//         Deltat_ctof_P_after_mmsq_cut->Fill(momentum, dt_ctof_p);
//         Deltat_ctof_pip_after_mmsq_cut->Fill(momentum, dt_ctof_pip);
// }
//
// void Histogram::Fill_deltat_after_mmsq_cut_Xpip(int pid, int charge, float dt_p, float dt_pip, float momentum) {
//         Deltat_P_after_mmsq_cut_singlepip->Fill(momentum, dt_p);
//         Deltat_pip_after_mmsq_cut_singlepip->Fill(momentum, dt_pip);
// }
// void Histogram::Fill_deltat_ctof_after_mmsq_cut_Xpip(int pid, int charge, float dt_ctof_p, float dt_ctof_pip,
//                                                      float momentum) {
//         Deltat_ctof_P_after_mmsq_cut_singlepip->Fill(momentum, dt_ctof_p);
//         Deltat_ctof_pip_after_mmsq_cut_singlepip->Fill(momentum, dt_ctof_pip);
// }
// void Histogram::Fill_deltat_vertex(int pid, int charge, float dt, float momentum) {
//         delta_t_vertex[0]->Fill(momentum, dt);
//         if (pid == ELECTRON) {
//                 delta_t_vertex[1]->Fill(momentum, dt);
//         } else {
//                 delta_t_vertex[2]->Fill(momentum, dt);
//         }
// }
// void Histogram::Fill_deltat_elect(int pid, int charge, float dt, float momentum) {
//         if (charge == -1) {
//                 delta_t_hist[0][1][0] -> Fill(momentum, dt);
//
//                 if (pid == ELECTRON) {
//                         delta_t_hist[0][1][1]->Fill(momentum, dt);
//                 } else {
//                         delta_t_hist[0][1][2]->Fill(momentum, dt);
//                 }
//         } else if (charge == 1) {
//                 delta_t_hist[0][0][0] -> Fill(momentum, dt);
//
//                 if (pid == -ELECTRON) {
//                         delta_t_hist[0][0][1]->Fill(momentum, dt);
//                 } else {
//                         delta_t_hist[0][0][2]->Fill(momentum, dt);
//                 }
//         }
// }
//
// void Histogram::Fill_deltat_prot(int pid, int charge, float dt, float momentum) {
//         //  for (size_t i = 0; i < with_id_num; i++) {
//         if (charge == 1) {
//                 delta_t_hist[2][0][0] -> Fill(momentum, dt);
//
//                 if (pid == PROTON) {
//                         delta_t_hist[2][0][1]->Fill(momentum, dt);
//
//                 } else {
//                         delta_t_hist[2][0][2]->Fill(momentum, dt);
//                 }
//         } else if (charge == -1) {
//                 delta_t_hist[2][1][0] -> Fill(momentum, dt);
//
//                 if (pid == -PROTON) {
//                         delta_t_hist[2][1][1]->Fill(momentum, dt);
//                 } else {
//                         delta_t_hist[2][1][2]->Fill(momentum, dt);
//                 }
//         }
// }
// void Histogram::Fill_deltat_pion(int pid, int charge, float dt, float momentum) {
//         //  for (size_t c = 0; c < charge_num; c++) {
//         //  for (size_t i = 0; i < with_id_num; i++) {
//         if (charge == 1) {
//                 delta_t_hist[1][0][0] -> Fill(momentum, dt);
//                 if (pid == PIP) {
//                         delta_t_hist[1][0][1]->Fill(momentum, dt);
//                 } else {
//                         delta_t_hist[1][0][2]->Fill(momentum, dt);
//                 }
//         } else if (charge == -1) {
//                 delta_t_hist[1][1][0] -> Fill(momentum, dt);
//                 if (pid == PIM) {
//                         delta_t_hist[1][1][1]->Fill(momentum, dt);
//                 } else
//                         delta_t_hist[1][1][2]->Fill(momentum, dt);
//         }
//
//         //}
// }
// void Histogram::Fill_deltat_kp(int pid, int charge, float dt, float momentum) {
//         //  for (size_t c = 0; c < charge_num; c++) {
//         //  for (size_t i = 0; i < with_id_num; i++) {
//         if (charge == 1) {
//                 delta_t_hist[3][0][0] -> Fill(momentum, dt);
//                 if (pid == KP) {
//                         delta_t_hist[3][0][1]->Fill(momentum, dt);
//                 } else {
//                         delta_t_hist[3][0][2]->Fill(momentum, dt);
//                 }
//         } else if (charge == -1) {
//                 delta_t_hist[3][1][0] -> Fill(momentum, dt);
//
//                 if (pid == KM) {
//                         delta_t_hist[3][1][1]->Fill(momentum, dt);
//                 } else
//                         delta_t_hist[3][1][2]->Fill(momentum, dt);
//         }
//
//         //}
// }
// void ::Histogram::Fill_ctof_e_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[0][1][1] -> Fill(momentum, dt_ctof); // e, -ve, without_cut
// }
// void ::Histogram::Fill_ctof_P_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[2][0][1] -> Fill(momentum, dt_ctof);
// }
// void ::Histogram::Fill_ctof_pip_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[1][0][1] -> Fill(momentum, dt_ctof);
// }
// void ::Histogram::Fill_ctof_pim_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[1][1][1] -> Fill(momentum, dt_ctof);
// }
// void ::Histogram::Fill_ctof_kp_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[3][0][1] -> Fill(momentum, dt_ctof);
// }
// void ::Histogram::Fill_ctof_km_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[3][1][1] -> Fill(momentum, dt_ctof);
// }
//
// void ::Histogram::Fill_ctof_e_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[0][1][0] -> Fill(momentum, dt_ctof);
// }
// void ::Histogram::Fill_ctof_P_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[2][0][0] -> Fill(momentum, dt_ctof);
// }
// void ::Histogram::Fill_ctof_pip_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[1][0][0] -> Fill(momentum, dt_ctof);
// }
// void ::Histogram::Fill_ctof_pim_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[1][1][0] -> Fill(momentum, dt_ctof);
// }
// void ::Histogram::Fill_ctof_kp_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[3][0][0] -> Fill(momentum, dt_ctof);
// }
// void ::Histogram::Fill_ctof_km_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
//         delta_t_hist_ctof[3][1][0] -> Fill(momentum, dt_ctof);
// }
// void Histogram::Write_deltat() {
//         for (size_t i = 0; i < with_id_num; i ++) {
//                 delta_t_vertex[i] -> SetXTitle("Momentum (GeV)");
//                 delta_t_vertex[i] -> SetYTitle("#Deltat");
//                 delta_t_vertex[i] -> SetOption("COLZ");
//                 delta_t_vertex[i] -> Write();
//                 delete delta_t_vertex[i];
//         }
//
//         for (size_t p = 0; p < particle_num; p ++) {
//                 for (size_t c = 0; c < charge_num; c ++) {
//                         for (size_t i = 0; i < with_id_num; i ++) {
//                                 delta_t_hist[p][c][i] -> SetXTitle("Momentum (GeV)");
//                                 delta_t_hist[p][c][i] -> SetYTitle("#Deltat");
//                                 delta_t_hist[p][c][i] -> SetOption("COLZ");
//                                 delta_t_hist[p][c][i] -> Write();
//                                 delete delta_t_hist[p][c][i];
//                         }
//                 }
//         }
//         for (size_t p = 0; p < particle_num; p ++) {
//                 for (size_t c = 0; c < charge_num; c ++) {
//                         for (size_t ct = 0; ct < cut_without_cut_num; ct ++) {
//                                 delta_t_hist_ctof[p][c][ct] -> SetXTitle("Momentum (GeV)");
//                                 delta_t_hist_ctof[p][c][ct] -> SetYTitle("#Deltat_ctof");
//                                 delta_t_hist_ctof[p][c][ct] -> SetOption("COLZ");
//                                 delta_t_hist_ctof[p][c][ct] -> Write();
//                                 delete delta_t_hist_ctof[p][c][ct];
//                         }
//                 }
//         }
//
//         Deltat_P_after_mmsq_cut->SetXTitle("Deltat_p_after_mmsq_cut");
//         Deltat_P_after_mmsq_cut->SetYTitle("Deltat_p_after_mmsq_cut");
//         Deltat_P_after_mmsq_cut->SetOption("COLZ");
//         Deltat_P_after_mmsq_cut->Write();
//         Deltat_pip_after_mmsq_cut->SetXTitle("Deltat_pip_after_mmsq_cut");
//         Deltat_pip_after_mmsq_cut->SetYTitle("Deltat_pip_after_mmsq_cut");
//         Deltat_pip_after_mmsq_cut->SetOption("COLZ");
//         Deltat_pip_after_mmsq_cut->Write();
//
//         Deltat_ctof_P_after_mmsq_cut->SetXTitle("Deltat_ctof_p_after_mmsq_cut");
//         Deltat_ctof_P_after_mmsq_cut->SetYTitle("Deltat_ctof_p_after_mmsq_cut");
//         Deltat_ctof_P_after_mmsq_cut->SetOption("COLZ");
//         Deltat_ctof_P_after_mmsq_cut->Write();
//         Deltat_ctof_pip_after_mmsq_cut->SetXTitle("Deltat_ctof_pip_after_mmsq_cut");
//         Deltat_ctof_pip_after_mmsq_cut->SetYTitle("Deltat_ctof_pip_after_mmsq_cut");
//         Deltat_ctof_pip_after_mmsq_cut->SetOption("COLZ");
//         Deltat_ctof_pip_after_mmsq_cut->Write();
//
//         Deltat_P_after_mmsq_cut_singlepip->SetXTitle("Deltat_p_after_mmsq_cut_singlepip");
//         Deltat_P_after_mmsq_cut_singlepip->SetYTitle("Deltat_p_after_mmsq_cut_singlepip");
//         Deltat_P_after_mmsq_cut_singlepip->SetOption("COLZ");
//         Deltat_P_after_mmsq_cut_singlepip->Write();
//         Deltat_pip_after_mmsq_cut_singlepip->SetXTitle("Deltat_pip_after_mmsq_cut_singlepip");
//         Deltat_pip_after_mmsq_cut_singlepip->SetYTitle("Deltat_pip_after_mmsq_cut_singlepip");
//         Deltat_pip_after_mmsq_cut_singlepip->SetOption("COLZ");
//         Deltat_pip_after_mmsq_cut_singlepip->Write();
//
//         Deltat_ctof_P_after_mmsq_cut_singlepip->SetXTitle("Deltat_ctof_p_after_mmsq_cut_singlepip");
//         Deltat_ctof_P_after_mmsq_cut_singlepip->SetYTitle("Deltat_ctof_p_after_mmsq_cut_singlepip");
//         Deltat_ctof_P_after_mmsq_cut_singlepip->SetOption("COLZ");
//         Deltat_ctof_P_after_mmsq_cut_singlepip->Write();
//         Deltat_ctof_pip_after_mmsq_cut_singlepip->SetXTitle("Deltat_ctof_pip_after_mmsq_cut_singlepip");
//         Deltat_ctof_pip_after_mmsq_cut_singlepip->SetYTitle("Deltat_ctof_pip_after_mmsq_cut_singlepip");
//         Deltat_ctof_pip_after_mmsq_cut_singlepip->SetOption("COLZ");
//         Deltat_ctof_pip_after_mmsq_cut_singlepip->Write();
// }
// // function below is for name and title of histogram vertex
// void Histogram::makeHists_MomVsBeta() {
//         for (size_t i = 0; i < with_id_num; i ++) {
//                 hname.append("mom_vs_beta_vertex");
//                 htitle.append("Momentum vs #beta vertex");
//                 hname.append("_");
//                 htitle.append(" ");
//                 hname.append(id_name[i]);
//                 htitle.append(id_name[i]);
//                 momvsbeta_vertex[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins,
//                 zero, 1.2); hname.clear(); htitle.clear();
//         }
//         //   // particle number = 4, 0e 1pi 2P 3 K
//         for (size_t p = 0; p < particle_num; p ++) {
//                 for (size_t c = 0; c < charge_num; c ++) {
//                         for (size_t i = 0; i < with_id_num; i ++) {
//                                 hname.append("mom_vs_beta_");
//                                 htitle.append("Momentum vs #beta ");
//                                 hname.append(particle_name[p]);
//                                 htitle.append(particle_name[p]);
//                                 hname.append("_");
//                                 htitle.append(" ");
//                                 hname.append(charge_name[c]);
//                                 htitle.append(charge_name[c]);
//                                 hname.append("_");
//                                 htitle.append(" ");
//                                 hname.append(id_name[i]);
//                                 htitle.append(id_name[i]);
//                                 momvsbeta_hist[p][c][i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min,
//                                 p_max, bins, zero, 1.2); hname.clear(); htitle.clear();
//                         }
//                 }
//         }
// }
// void Histogram::Fill_MomVsBeta_vertex(int pid, int charge, float P, float beta) {
//         if (beta != 0) {
//                 momvsbeta_vertex[0] -> Fill(P, beta);
//                 if (pid == ELECTRON) {
//                         momvsbeta_vertex[1] -> Fill(P, beta);
//
//                 } else {
//                         momvsbeta_vertex[2] -> Fill(P, beta);
//                 }
//         }
// }
//
// void Histogram::Fill_MomVsBeta(int pid, int charge, float P, float beta) {
//         int good_ID = 0;
//         if (beta != 0) {
//                 momentum->Fill(P);
//                 for (size_t p = 0; p < particle_num; p ++) {
//                         switch (p) {
//                         case 0:
//                                 good_ID = -ELECTRON;
//                                 break;
//                         case 1:
//                                 good_ID = PIP;
//                                 break;
//                         case 2:
//                                 good_ID = PROTON;
//                                 break;
//                         case 3:
//                                 good_ID = KP;
//                                 break;
//                         }
//
//                         /*momvsbeta_hist[p][0][0]->Fill(P, beta);
//                            if (good_ID == abs(pid)) {
//                            momvsbeta_hist[p][0][1]->Fill(P, beta);
//                            } else {
//                            momvsbeta_hist[p][0][2]->Fill(P, beta);
//                            }*/
//
//                         if (charge == -1) {
//                                 momvsbeta_hist[p][1][0] -> Fill(P, beta);
//                                 if (-good_ID == pid) { // - good_ID thyo paila
//                                         momvsbeta_hist[p][1][1]->Fill(P, beta);
//                                 } else {
//                                         momvsbeta_hist[p][1][2]->Fill(P, beta);
//                                 }
//                         } else if (charge == 1) {
//                                 momvsbeta_hist[p][0][0] -> Fill(P, beta);
//                                 if (good_ID == pid) {
//                                         momvsbeta_hist[p][0][1]->Fill(P, beta);
//                                 } else {
//                                         momvsbeta_hist[p][0][2]->Fill(P, beta);
//                                 }
//                         }
//                 }
//         }
// }
//
// void Histogram::Write_MomVsBeta() {
//         for (size_t i = 0; i < with_id_num; i++) {
//                 momvsbeta_vertex[i]->SetXTitle("Momentum (GeV)");
//                 momvsbeta_vertex[i]->SetYTitle("#beta");
//                 momvsbeta_vertex[i]->SetOption("COLZ");
//                 momvsbeta_vertex[i]->Write();
//                 delete momvsbeta_vertex[i];
//         }
//
//         momentum->SetXTitle("Momentum (GeV)");
//         momentum->Write();
//         // delete momentum;
//         for (size_t p = 0; p < particle_num; p++) {
//                 for (size_t c = 0; c < charge_num; c++) {
//                         for (size_t i = 0; i < with_id_num; i++) {
//                                 momvsbeta_hist[p][c][i]->SetXTitle("Momentum (GeV)");
//                                 momvsbeta_hist[p][c][i]->SetYTitle("#beta");
//                                 momvsbeta_hist[p][c][i]->SetOption("COLZ");
//                                 momvsbeta_hist[p][c][i]->Write();
//                                 delete momvsbeta_hist[p][c][i];
//                         }
//                 }
//         }
// }
void Histogram::Write_EC() {
  for (int i = 0; i < sec_num; i++) {
    EC_sampling_fraction[i]->SetXTitle("Momentum (GeV)");
    EC_sampling_fraction[i]->SetYTitle("Sampling Fraction");
    EC_sampling_fraction[i]->SetOption("COLZ");
    EC_sampling_fraction[i]->Write();
    delete EC_sampling_fraction[i];

    PCAL_VS_ECAL[i]->SetXTitle("Pcal (GeV)");
    PCAL_VS_ECAL[i]->SetYTitle("Ecal_(Ecin+Ecout)");
    PCAL_VS_ECAL[i]->SetOption("COLZ");
    PCAL_VS_ECAL[i]->Write();
    delete PCAL_VS_ECAL[i];

    for (int n = 0; n < slice_num; n++) {
      sf_[i][n]->SetXTitle("Momentum (GeV)");
      sf_[i][n]->Write();
      delete sf_[i][n];
    }
  }

  // lu_side_distribution->SetXTitle("side_U");
  // lu_side_distribution->Write();
  // lv_side_distribution->SetXTitle("side_V");
  // lv_side_distribution->Write();
  // lw_side_distribution->SetXTitle("side_W");
  // lw_side_distribution->Write();
  for (int j = 0; j < cut_y_n; j++) {
    PCAL_FID_CUT[j]->SetXTitle("x/cm");
    PCAL_FID_CUT[j]->SetYTitle("y/cm");
    PCAL_FID_CUT[j]->SetOption("COLZ");
    PCAL_FID_CUT[j]->Write();
    delete PCAL_FID_CUT[j];
    DCr1_FID_CUT[j]->SetXTitle("x/cm");
    DCr1_FID_CUT[j]->SetYTitle("y/cm");
    DCr1_FID_CUT[j]->SetOption("COLZ");
    DCr1_FID_CUT[j]->Write();
    delete DCr1_FID_CUT[j];
    DCr2_FID_CUT[j]->SetXTitle("x/cm");
    DCr2_FID_CUT[j]->SetYTitle("y/cm");
    DCr2_FID_CUT[j]->SetOption("COLZ");
    DCr2_FID_CUT[j]->Write();
    delete DCr2_FID_CUT[j];
    DCr3_FID_CUT[j]->SetXTitle("x/cm");
    DCr3_FID_CUT[j]->SetYTitle("y/cm");
    DCr3_FID_CUT[j]->SetOption("COLZ");
    DCr3_FID_CUT[j]->Write();
    delete DCr3_FID_CUT[j];
  }
  vertex_vz->SetXTitle("vertex_vz");
  vertex_vz->Write();
  // delete vertex_vz;/
  // pid_size1->SetXTitle("pid_size");
  // pid_size1->Write();
  // pid_check->SetXTitle("pid");
  // pid_check->Write();
  // mom_diff_e_pim->SetXTitle("mom_diff_e_pim");
  // mom_diff_e_pim->Write();
  // mom_diff_p_pip->SetXTitle("mom_diff_p_pip");
  // mom_diff_p_pip->Write();
  // mom_diff_e_pim_x->SetXTitle("mom_diff_e_pim_x");
  // mom_diff_e_pim_x->Write();
  // mom_diff_p_pip_x->SetXTitle("mom_diff_p_pip_x");
  // mom_diff_p_pip_x->Write();
  // mom_diff_e_pim_y->SetXTitle("mom_diff_e_pim_y");
  // mom_diff_e_pim_y->Write();
  // mom_diff_p_pip_y->SetXTitle("mom_diff_p_pip_y");
  // mom_diff_p_pip_y->Write();
  // mom_diff_e_pim_z->SetXTitle("mom_diff_e_pim_z");
  // mom_diff_e_pim_z->Write();
  // mom_diff_p_pip_z->SetXTitle("mom_diff_p_pip_z");
  // mom_diff_p_pip_z->Write();
  //
  // mom_e->SetXTitle("mom_e");
  // mom_e->Write();
  // mom_p->SetXTitle("mom_p");
  // mom_p->Write();
  // mom_e_z->SetXTitle("mom_e_z");
  // mom_e_z->Write();
  // mom_p_z->SetXTitle("mom_p_z");
  // mom_p_z->Write();
  // mom_pip->SetXTitle("mom_pip");
  // mom_pip->Write();
  // mom_pim->SetXTitle("mom_pim");
  // mom_pim->Write();
  // mom_pip_z->SetXTitle("mom_pip_z");
  // mom_pip_z->Write();
  // mom_pim_z->SetXTitle("mom_pim_z");
  // mom_pim_z->Write();
}
