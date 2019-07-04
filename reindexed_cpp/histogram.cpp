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
  momentum = new TH1D("mom", "mom", bins, p_min, p_max);

  vertex_vz = new TH1D("vertex_position", "vertex_position", bins, -40, 40);

  theta_P_vs_mass_pip_pim =
      new TH2D("theta_p_vs_inv_mass_pip_pim", "theta_P_vs_inv_mass_pip_pim", bins, 0.2, 2.2, 100, 0, 180);
  theta_pim_vs_mass_Ppip =
      new TH2D("theta_pim_vs_inv_mass_P_pip", "theta_pim_vs_inv_mass_P_pip", bins, 1, 3.5, 100, 0, 180);
  theta_pip_vs_mass_Ppim =
      new TH2D("theta_pip_vs_inv_mass_P_pim", "theta_pip_vs_inv_mass_P_pim", bins, 1, 3.5, 100, 0, 180);

  theta_P_lab_vs_mass_pip_pim =
      new TH2D("theta_P_lab_vs_inv_mass_pip_pim", "theta_P_lab_vs_inv_mass_pip_pim", bins, 0.2, 2.2, 100, 0, 75);
  theta_pim_lab_vs_mass_Ppip =
      new TH2D("theta_pim_lab_vs_inv_mass_P_pip", "theta_pim_lab_vs_inv_mass_P_pip", bins, 1, 3.5, 100, 0, 125);
  theta_pip_lab_vs_mass_Ppim =
      new TH2D("theta_pip_lab_vs_inv_mass_P_pim", "theta_pip_lab_vs_inv_mass_P_pim", bins, 1, 3.5, 100, 0, 125);
  lu_side_distribution = new TH1D("lu_side_distribution", "lu_side_distribution", 50, 0, 400);
  lv_side_distribution = new TH1D("lv_side_distribution", "lv_side_distribution", 50, 0, 450);
  lw_side_distribution = new TH1D("lw_side_distribution", "lw_side_distribution", 50, 0, 450);

  W_hist_twopi_all_sec = new TH1D("W_twoPions_all_sec", "W_twoPions_all_sec", bins, 1, 3.5);
  inv_mass_P_pip_all_sec = new TH1D("invariant_mass_P_pip_all_sec", "invariant_mass_P_pip_all_sec", bins, 1, 2.8);
  inv_mass_P_pim_all_sec = new TH1D("invariant_mass_P_pim_all_sec", "invariant_mass_P_pim_all_sec", bins, 1, 2.6);
  inv_mass_pip_pim_all_sec =
      new TH1D("invariant_mass_pip_pim_all_sec", "invariant_mass_pip_pim_all_sec", bins, 0.2, 1.8);
  W_hist_Xpip_all_sec = new TH1D("W_Xpip_all_sec", "W_Xpip_all_sec", bins, 0.9, 2.8);

  makeHists_deltat();
  // makeHists_MomVsBeta();
  makeHists_WvsQ2();
  makeHists_MM();
  // Make_hist_cc();
  makeHists_EC_sf();
  makeHists_pcal_fid_cuts();
}

Histogram::~Histogram() {}
// W and Q^2

float Histogram::mm_lim_min(int mm_number, int mm_events_number) {
  if (mm_number == 0 && mm_events_number < 4) {
    return -2.50;
  } else if (mm_number == 0 && mm_events_number > 4) {
    return -2.0;
  } else if (mm_number == 0 && mm_events_number == 4) {
    return 0;
  } else if (mm_number == 1 && mm_events_number < 4 && mm_events_number != 1) {
    return -1;
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
    return 4.50;
  } else if (mm_number == 0 && mm_events_number > 4) {
    return 6.50;
  } else if (mm_number == 0 && mm_events_number == 4) {
    return 6.0;
  } else if (mm_number == 1 && mm_events_number < 4 && mm_events_number != 1) {
    return 3.7;
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
    PCAL_VS_ECAL[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, zero, 1, bins, zero, 1);
    hname.clear();
    htitle.clear();
  }
}
void Histogram::makeHists_WvsQ2() {
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

void Histogram::Fill_WvsQ2(double W, double Q2, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_vs_Q2[sec_number]->Fill(W, Q2);
    W_hist[sec_number]->Fill(W);

    Q2_hist[sec_number]->Fill(Q2);
  }
}
void Histogram::Fill_EC_sampling_fraction(double momentum, double sf, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    EC_sampling_fraction[sec_number]->Fill(momentum, sf);
  }
}

void Histogram::Fill_PCAL_VS_ECAL(float pcal, float ecal, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    PCAL_VS_ECAL[sec_number]->Fill(pcal, ecal);
  }
}
void Histogram::Fill_WvsmmSQ_ep(double W, double mmSQ, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_hist_ep[sec_number]->Fill(W);
    W_vs_mmSQ_ep[sec_number][0]->Fill(W, mmSQ);
  }
}
void Histogram::Fill_W_2pi_all_sec(double W, double W_dpp, double delta_zero_, double rho_) {
  W_hist_twopi_all_sec->Fill(W);
  inv_mass_P_pip_all_sec->Fill(W_dpp);
  inv_mass_P_pim_all_sec->Fill(delta_zero_);
  inv_mass_pip_pim_all_sec->Fill(rho_);
}
void Histogram::Fill_W_hist_Xpip_all_sec(double W) { W_hist_Xpip_all_sec->Fill(W); }
void Histogram::Fill_WvsmmSQ_2pi(double W, double W_dpp, double delta_zero_, double rho_, double mmSQ, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_hist_2pi[sec_number]->Fill(W);
    W_hist_delta_pp[sec_number]->Fill(W_dpp);
    W_hist_delta_zero[sec_number]->Fill(delta_zero_);
    W_hist_rho[sec_number]->Fill(rho_);
    W_vs_mmSQ_2pi[sec_number][0]->Fill(W, mmSQ);
  }
}
void Histogram::Fill_WvsmmSQ_singlepip(double W, double mmSQ, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_hist_singlepip[sec_number]->Fill(W);
    W_vs_mmSQ_singlepip[sec_number][0]->Fill(W, mmSQ);
  }
}
void Histogram::Fill_WvsmmSQ_anti_ep(double W, double mmSQ, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_vs_mmSQ_ep[sec_number][1]->Fill(W, mmSQ);
  }
}
void Histogram::Fill_WvsmmSQ_anti_2pi(double W, double W_dpp, double delta_zero_, double rho_, double mmSQ,
                                      int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_vs_mmSQ_2pi[sec_number][1]->Fill(W, mmSQ);
  }
}
void Histogram::Fill_WvsmmSQ_anti_singlepip(double W, double mmSQ, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_vs_mmSQ_singlepip[sec_number][1]->Fill(W, mmSQ);
  }
}

// void Histogram::Make_hist_cc() {
//   for (int i = 0; i < cc_num; i++) {
//     hname.append("cc_total_");
//     htitle.append("cc_total_");
//     hname.append(cc_name[i]);
//     htitle.append(cc_name[i]);
//     cherenkov_total[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 60);
//     hname.clear();
//     htitle.clear();
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
// void Histogram::Fill_hist_cc_tot(float tot_el) { cherenkov_total[0]->Fill(tot_el); }
// void Histogram::Fill_hist_cc_ltcc(float ltcc_el) { cherenkov_ltcc[0]->Fill(ltcc_el); }
// void Histogram::Fill_hist_cc_htcc(float htcc_el) { cherenkov_htcc[0]->Fill(htcc_el); }
// void Histogram::Fill_hist_cc_tot_pim(float tot_pim) { cherenkov_total[1]->Fill(tot_pim); }
// void Histogram::Fill_hist_cc_ltcc_pim(float ltcc_pim) { cherenkov_ltcc[1]->Fill(ltcc_pim); }
// void Histogram::Fill_hist_cc_htcc_pim(float htcc_pim) { cherenkov_htcc[1]->Fill(htcc_pim); }
// void Histogram::Fill_hist_cc_tot_pip(float tot_pip) { cherenkov_total[2]->Fill(tot_pip); }
// void Histogram::Fill_hist_cc_ltcc_pip(float ltcc_pip) { cherenkov_ltcc[2]->Fill(ltcc_pip); }
// void Histogram::Fill_hist_cc_htcc_pip(float htcc_pip) { cherenkov_htcc[2]->Fill(htcc_pip); }
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

void Histogram::Fill_theta_P_inv_mass(float inv_mass, float theta) { theta_P_vs_mass_pip_pim->Fill(inv_mass, theta); }
void Histogram::Fill_theta_pim_inv_mass(float inv_mass, float theta) { theta_pim_vs_mass_Ppip->Fill(inv_mass, theta); }
void Histogram::Fill_theta_pip_inv_mass(float inv_mass, float theta) { theta_pip_vs_mass_Ppim->Fill(inv_mass, theta); }
void Histogram::Fill_theta_P_lab_inv_mass(float inv_mass, float theta_lab) {
  theta_P_lab_vs_mass_pip_pim->Fill(inv_mass, theta_lab);
}
void Histogram::Fill_theta_pim_lab_inv_mass(float inv_mass, float theta_lab) {
  theta_pim_lab_vs_mass_Ppip->Fill(inv_mass, theta_lab);
}
void Histogram::Fill_theta_pip_lab_inv_mass(float inv_mass, float theta_lab) {
  theta_pip_lab_vs_mass_Ppim->Fill(inv_mass, theta_lab);
}
void Histogram::Fill_lu_dist(float li) { lu_side_distribution->Fill(li); }
void Histogram::Fill_lv_dist(float li) { lv_side_distribution->Fill(li); }
void Histogram::Fill_lw_dist(float li) { lw_side_distribution->Fill(li); }

void Histogram::Fill_vertex_vz(float vz) { vertex_vz->Fill(vz); }
// void Histogram::Fill_theta_P(float theta_p, float theta_pip_, float theta_pim_) {
//         theta_prot->Fill(theta_p);
//         theta_pip->Fill(theta_pip_);
//         theta_pim->Fill(theta_pim_);
// }
// void Histogram::Fill_Phi_cm(float Phi_p, float Phi_pip_, float Phi_pim_) {
//         Phi_prot->Fill(Phi_p);
//         Phi_pip->Fill(Phi_pip_);
//         Phi_pim->Fill(Phi_pim_);
// }

void Histogram::Write_WvsQ2() {
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
  W_hist_twopi_all_sec->SetXTitle("W_P_pip_pim (GeV)");
  W_hist_twopi_all_sec->Write();
  // delete W_hist_twopi_all_sec;
  inv_mass_P_pip_all_sec->SetXTitle("inv_mass_P_pip (GeV)");
  inv_mass_P_pip_all_sec->Write();
  // delete inv_mass_P_pip_all_sec;
  inv_mass_P_pim_all_sec->SetXTitle("inv_mass_P_pim (GeV)");
  inv_mass_P_pim_all_sec->Write();
  // delete inv_mass_P_pip_all_sec;
  inv_mass_pip_pim_all_sec->SetXTitle("inv_mass_pip_pim (GeV)");
  inv_mass_pip_pim_all_sec->Write();
  // delete invariant_mass_pip_pim_all_sec;
  W_hist_Xpip_all_sec->SetXTitle("W_Xpip_all_sec (GeV)");
  W_hist_Xpip_all_sec->Write();
  // delete W_hist_Xpip_all_sec;
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

void Histogram::Fill_ep_mm(double mm, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][0][sec_number]->Fill(mm);
  }
}
void Histogram::Fill_ep_mmSQ(double mm, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][0][sec_number]->Fill(mm);
  }
}
void Histogram::Fill_2pion_mm(double mm, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][1][sec_number]->Fill(mm);
  }
}
void Histogram::Fill_2pion_mmSQ(double mm, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][1][sec_number]->Fill(mm);
  }
}
void Histogram::Fill_pip_mm(double mm, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][2][sec_number]->Fill(mm);
  }
}
void Histogram::Fill_pip_mmSQ(double mm, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][2][sec_number]->Fill(mm);
  }
}
void Histogram::Fill_pim_mm(double mm, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][3][sec_number]->Fill(mm);
  }
}
void Histogram::Fill_pim_mmSQ(double mm, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][3][sec_number]->Fill(mm);
  }
}
void Histogram::Fill_MM_wop_e_prime(double mm_1, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][4][sec_number]->Fill(mm_1);
  }
}
void Histogram::Fill_MMSQ_wop_e_prime(double mm_1, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][4][sec_number]->Fill(mm_1);
  }
}
void Histogram::Fill_MM_wop_2pion(double mm_1, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][5][sec_number]->Fill(mm_1);
  }
}
void Histogram::Fill_MMSQ_wop_2pion(double mm_1, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][5][sec_number]->Fill(mm_1);
  }
}
void Histogram::Fill_MM_wop_pip(double mm_1, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][6][sec_number]->Fill(mm_1);
  }
}
void Histogram::Fill_MMSQ_wop_pip(double mm_1, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][6][sec_number]->Fill(mm_1);
  }
}
void Histogram::Fill_MM_wop_pim(double mm_1, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[0][7][sec_number]->Fill(mm_1);
  }
}
void Histogram::Fill_MMSQ_wop_pim(double mm_1, int sec_number) {
  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    MM_hist[1][7][sec_number]->Fill(mm_1);
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
void Histogram::makeHists_deltat() {
  for (size_t i = 0; i < with_id_num; i++) {
    hname.append("delta_t_vertex");
    htitle.append("#Deltat vertex particle");
    hname.append("_");
    htitle.append(" ");
    hname.append(id_name[i]);
    htitle.append(id_name[i]);
    delta_t_vertex[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, Dt_min, Dt_max);
    hname.clear();
    htitle.clear();
  }

  for (size_t p = 0; p < particle_num; p++) {
    for (size_t c = 0; c < charge_num; c++) {
      for (size_t i = 0; i < with_id_num; i++) {
        hname.append("delta_t_");
        htitle.append("#Deltat ");
        hname.append(particle_name[p]);
        htitle.append(particle_name[p]);
        hname.append("_");
        htitle.append(" ");
        hname.append(charge_name[c]);
        htitle.append(charge_name[c]);
        hname.append("_");
        htitle.append(" ");
        hname.append(id_name[i]);
        htitle.append(id_name[i]);
        delta_t_hist[p][c][i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, Dt_min, Dt_max);
        hname.clear();
        htitle.clear();
      }
    }
  }
  for (size_t p = 0; p < particle_num; p++) {
    for (size_t c = 0; c < charge_num; c++) {
      for (size_t ct = 0; ct < cut_without_cut_num; ct++) {
        hname.append("delta_t_");
        htitle.append("#Deltat ");
        hname.append(particle_name[p]);
        htitle.append(particle_name[p]);
        hname.append("_ctof_");
        htitle.append(" ctof ");
        hname.append(charge_name[c]);
        htitle.append(charge_name[c]);
        hname.append("_");
        htitle.append(" ");
        hname.append(cut_without_cut_name[ct]);
        htitle.append(cut_without_cut_name[ct]);
        delta_t_hist_ctof[p][c][ct] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, Dt_min, Dt_max);
        hname.clear();
        htitle.clear();
      }
    }
  }
}

void Histogram::Fill_deltat_vertex(int pid, int charge, float dt, float momentum) {
  delta_t_vertex[0]->Fill(momentum, dt);
  if (pid == ELECTRON) {
    delta_t_vertex[1]->Fill(momentum, dt);
  } else {
    delta_t_vertex[2]->Fill(momentum, dt);
  }
}
void Histogram::Fill_deltat_elect(int pid, int charge, float dt, float momentum) {
  if (charge == -1) {
    delta_t_hist[0][1][0]->Fill(momentum, dt);

    if (pid == ELECTRON) {
      delta_t_hist[0][1][1]->Fill(momentum, dt);
    } else {
      delta_t_hist[0][1][2]->Fill(momentum, dt);
    }
  } else if (charge == 1) {
    delta_t_hist[0][0][0]->Fill(momentum, dt);

    if (pid == -ELECTRON) {
      delta_t_hist[0][0][1]->Fill(momentum, dt);
    } else {
      delta_t_hist[0][0][2]->Fill(momentum, dt);
    }
  }
}

void Histogram::Fill_deltat_prot(int pid, int charge, float dt, float momentum) {
  //  for (size_t i = 0; i < with_id_num; i++) {
  if (charge == 1) {
    delta_t_hist[2][0][0]->Fill(momentum, dt);

    if (pid == PROTON) {
      delta_t_hist[2][0][1]->Fill(momentum, dt);

    } else {
      delta_t_hist[2][0][2]->Fill(momentum, dt);
    }
  } else if (charge == -1) {
    delta_t_hist[2][1][0]->Fill(momentum, dt);

    if (pid == -PROTON) {
      delta_t_hist[2][1][1]->Fill(momentum, dt);
    } else {
      delta_t_hist[2][1][2]->Fill(momentum, dt);
    }
  }
}
void Histogram::Fill_deltat_pip(int pid, int charge, float dt, float momentum) {
  //  for (size_t c = 0; c < charge_num; c++) {
  //  for (size_t i = 0; i < with_id_num; i++) {
  if (charge == 1) {
    delta_t_hist[1][0][0]->Fill(momentum, dt);
    if (pid == PIP) {
      delta_t_hist[1][0][1]->Fill(momentum, dt);
    } else {
      delta_t_hist[1][0][2]->Fill(momentum, dt);
    }
  } else if (charge == -1) {
    delta_t_hist[1][1][0]->Fill(momentum, dt);
    if (pid == PIM) {
      delta_t_hist[1][1][1]->Fill(momentum, dt);
    } else
      delta_t_hist[1][1][2]->Fill(momentum, dt);
  }

  //}
}
void Histogram::Fill_deltat_kp(int pid, int charge, float dt, float momentum) {
  //  for (size_t c = 0; c < charge_num; c++) {
  //  for (size_t i = 0; i < with_id_num; i++) {
  if (charge == 1) {
    delta_t_hist[3][0][0]->Fill(momentum, dt);
    if (pid == KP) {
      delta_t_hist[3][0][1]->Fill(momentum, dt);
    } else {
      delta_t_hist[3][0][2]->Fill(momentum, dt);
    }
  } else if (charge == -1) {
    delta_t_hist[3][1][0]->Fill(momentum, dt);

    if (pid == KM) {
      delta_t_hist[3][1][1]->Fill(momentum, dt);
    } else
      delta_t_hist[3][1][2]->Fill(momentum, dt);
  }

  //}
}
void ::Histogram::Fill_ctof_e_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[0][1][1]->Fill(momentum, dt_ctof);  // e, -ve, without_cut
}
void ::Histogram::Fill_ctof_P_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[2][0][1]->Fill(momentum, dt_ctof);
}
void ::Histogram::Fill_ctof_pip_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[1][0][1]->Fill(momentum, dt_ctof);
}
void ::Histogram::Fill_ctof_pim_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[1][1][1]->Fill(momentum, dt_ctof);
}
void ::Histogram::Fill_ctof_kp_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[3][0][1]->Fill(momentum, dt_ctof);
}
void ::Histogram::Fill_ctof_km_without_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[3][1][1]->Fill(momentum, dt_ctof);
}

void ::Histogram::Fill_ctof_e_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[0][1][0]->Fill(momentum, dt_ctof);
}
void ::Histogram::Fill_ctof_P_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[2][0][0]->Fill(momentum, dt_ctof);
}
void ::Histogram::Fill_ctof_pip_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[1][0][0]->Fill(momentum, dt_ctof);
}
void ::Histogram::Fill_ctof_pim_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[1][1][0]->Fill(momentum, dt_ctof);
}
void ::Histogram::Fill_ctof_kp_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[3][0][0]->Fill(momentum, dt_ctof);
}
void ::Histogram::Fill_ctof_km_with_cut_hist(int pid, int charge, float dt_ctof, float momentum) {
  delta_t_hist_ctof[3][1][0]->Fill(momentum, dt_ctof);
}
void Histogram::Write_deltat() {
  for (size_t i = 0; i < with_id_num; i++) {
    delta_t_vertex[i]->SetXTitle("Momentum (GeV)");
    delta_t_vertex[i]->SetYTitle("#Deltat");
    delta_t_vertex[i]->SetOption("COLZ");
    delta_t_vertex[i]->Write();
    delete delta_t_vertex[i];
  }

  for (size_t p = 0; p < particle_num; p++) {
    for (size_t c = 0; c < charge_num; c++) {
      for (size_t i = 0; i < with_id_num; i++) {
        delta_t_hist[p][c][i]->SetXTitle("Momentum (GeV)");
        delta_t_hist[p][c][i]->SetYTitle("#Deltat");
        delta_t_hist[p][c][i]->SetOption("COLZ");
        delta_t_hist[p][c][i]->Write();
        delete delta_t_hist[p][c][i];
      }
    }
  }
  for (size_t p = 0; p < particle_num; p++) {
    for (size_t c = 0; c < charge_num; c++) {
      for (size_t ct = 0; ct < cut_without_cut_num; ct++) {
        delta_t_hist_ctof[p][c][ct]->SetXTitle("Momentum (GeV)");
        delta_t_hist_ctof[p][c][ct]->SetYTitle("#Deltat_ctof");
        delta_t_hist_ctof[p][c][ct]->SetOption("COLZ");
        delta_t_hist_ctof[p][c][ct]->Write();
        delete delta_t_hist_ctof[p][c][ct];
      }
    }
  }
}
// function below is for name and title of histogram vertex
// void Histogram::makeHists_MomVsBeta() {
//   for (size_t i = 0; i < with_id_num; i++) {
//     hname.append("mom_vs_beta_vertex");
//     htitle.append("Momentum vs #beta vertex");
//     hname.append("_");
//     htitle.append(" ");
//     hname.append(id_name[i]);
//     htitle.append(id_name[i]);
//     momvsbeta_vertex[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, zero, 1.2);
//     hname.clear();
//     htitle.clear();
//   }
//   // particle number = 4, 0e 1pi 2P 3 K
//   for (size_t p = 0; p < particle_num; p++) {
//     for (size_t c = 0; c < charge_num; c++) {
//       for (size_t i = 0; i < with_id_num; i++) {
//         hname.append("mom_vs_beta_");
//         htitle.append("Momentum vs #beta ");
//         hname.append(particle_name[p]);
//         htitle.append(particle_name[p]);
//         hname.append("_");
//         htitle.append(" ");
//         hname.append(charge_name[c]);
//         htitle.append(charge_name[c]);
//         hname.append("_");
//         htitle.append(" ");
//         hname.append(id_name[i]);
//         htitle.append(id_name[i]);
//         momvsbeta_hist[p][c][i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, zero, 1.2);
//         hname.clear();
//         htitle.clear();
//       }
//     }
//   }
// }
// void Histogram::Fill_MomVsBeta_vertex(int pid, int charge, double P, double beta) {
//   if (beta != 0) {
//     momvsbeta_vertex[0]->Fill(P, beta);
//     if (pid == ELECTRON) {
//       momvsbeta_vertex[1]->Fill(P, beta);
//
//     } else {
//       momvsbeta_vertex[2]->Fill(P, beta);
//     }
//   }
// }

// void Histogram::Fill_MomVsBeta(int pid, int charge, double P, double beta) {
//   int good_ID = 0;
//   if (beta != 0) {
//     momentum->Fill(P);
//     for (size_t p = 0; p < particle_num; p++) {
//       switch (p) {
//         case 0:
//           good_ID = -ELECTRON;
//           break;
//         case 1:
//           good_ID = PIP;
//           break;
//         case 2:
//           good_ID = PROTON;
//           break;
//         case 3:
//           good_ID = KP;
//           break;
//       }
//
//       /*momvsbeta_hist[p][0][0]->Fill(P, beta);
//          if (good_ID == abs(pid)) {
//          momvsbeta_hist[p][0][1]->Fill(P, beta);
//          } else {
//          momvsbeta_hist[p][0][2]->Fill(P, beta);
//          }*/
//
//       if (charge == -1) {
//         momvsbeta_hist[p][1][0]->Fill(P, beta);
//         if (-good_ID == pid) {  // - good_ID thyo paila
//           momvsbeta_hist[p][1][1]->Fill(P, beta);
//         } else {
//           momvsbeta_hist[p][1][2]->Fill(P, beta);
//         }
//       } else if (charge == 1) {
//         momvsbeta_hist[p][0][0]->Fill(P, beta);
//         if (good_ID == pid) {
//           momvsbeta_hist[p][0][1]->Fill(P, beta);
//         } else {
//           momvsbeta_hist[p][0][2]->Fill(P, beta);
//         }
//       }
//     }
//   }
// }

// void Histogram::Write_MomVsBeta() {
//   for (size_t i = 0; i < with_id_num; i++) {
//     momvsbeta_vertex[i]->SetXTitle("Momentum (GeV)");
//     momvsbeta_vertex[i]->SetYTitle("#beta");
//     momvsbeta_vertex[i]->SetOption("COLZ");
//     momvsbeta_vertex[i]->Write();
//     delete momvsbeta_vertex[i];
//   }
//
//   momentum->SetXTitle("Momentum (GeV)");
//   momentum->Write();
//   // delete momentum;
//   for (size_t p = 0; p < particle_num; p++) {
//     for (size_t c = 0; c < charge_num; c++) {
//       for (size_t i = 0; i < with_id_num; i++) {
//         momvsbeta_hist[p][c][i]->SetXTitle("Momentum (GeV)");
//         momvsbeta_hist[p][c][i]->SetYTitle("#beta");
//         momvsbeta_hist[p][c][i]->SetOption("COLZ");
//         momvsbeta_hist[p][c][i]->Write();
//         delete momvsbeta_hist[p][c][i];
//       }
//     }
//   }
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
  // delete vertex_vz;

  theta_pim_vs_mass_Ppip->SetXTitle("Inv_m_P_pi+");
  theta_pim_vs_mass_Ppip->SetYTitle("theta_pi_cm");
  theta_pim_vs_mass_Ppip->SetOption("COLZ");
  theta_pim_vs_mass_Ppip->Write();
  // delete theta_pim_vs_mass_Ppip;

  theta_pip_vs_mass_Ppim->SetXTitle("Inv_m_P-pi-");
  theta_pip_vs_mass_Ppim->SetYTitle("theta_pi+_cm");
  theta_pip_vs_mass_Ppim->SetOption("COLZ");
  theta_pip_vs_mass_Ppim->Write();
  // delete theta_pip_vs_mass_Ppim;

  theta_P_vs_mass_pip_pim->SetXTitle("Inv_m_Pi+pi-");
  theta_P_vs_mass_pip_pim->SetYTitle("theta_P_cm");
  theta_P_vs_mass_pip_pim->SetOption("COLZ");
  theta_P_vs_mass_pip_pim->Write();
  // delete theta_P_vs_mass_pip_pim;

  theta_pim_lab_vs_mass_Ppip->SetXTitle("Inv_m_P_pi+");
  theta_pim_lab_vs_mass_Ppip->SetYTitle("theta_pi-_lab");
  theta_pim_lab_vs_mass_Ppip->SetOption("COLZ");
  theta_pim_lab_vs_mass_Ppip->Write();
  // delete theta_pim_lab_vs_mass_Ppip;

  theta_pip_lab_vs_mass_Ppim->SetXTitle("Inv_m_P-pi-");
  theta_pip_lab_vs_mass_Ppim->SetYTitle("theta_pi+_lab");
  theta_pip_lab_vs_mass_Ppim->SetOption("COLZ");
  theta_pip_lab_vs_mass_Ppim->Write();
  // delete theta_pip_lab_vs_mass_Ppim;

  theta_P_lab_vs_mass_pip_pim->SetXTitle("Inv_m_Pi+pi-");
  theta_P_lab_vs_mass_pip_pim->SetYTitle("theta_P_lab");
  theta_P_lab_vs_mass_pip_pim->SetOption("COLZ");
  theta_P_lab_vs_mass_pip_pim->Write();
}
