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

  // MM_neutron =
  //     new TH1D("missMassNEutron", "missMassNeutron", bins, zero, p_max);
  theta_prot = new TH1D("theta_prot_dist_cm", "theta_dist-prot_cm", bins, zero, 180);
  theta_pip = new TH1D("theta_pip_dist_cm", "theta_dist-pip_cm", bins, zero, 180);
  theta_pim = new TH1D("theta_pim_dist_cm", "theta_dist-pim_cm", bins, zero, 180);
  Phi_prot = new TH1D("Phi_prot_dist_cm", "Phi_dist-prot_cm", bins, -180, 180);
  Phi_pip = new TH1D("Phi_pip_dist_cm", "Phi_dist-pip_cm", bins, -180, 180);
  Phi_pim = new TH1D("Phi_pim_dist_cm", "Phi_dist-pim_cm", bins, -180, 180);

  vertex_vz = new TH1D("vertex_position", "vertex_position", bins, -40, 40);

  theta_vs_phi_cm = new TH2D("theta_vs_phi_cm", "theta_vs_phi_cm", bins, zero, 60, 100, -180, 180);
  sf_vs_lv_distance_on_v_side = new TH2D("E/p_vs_lv", "E/p_vs_lv", bins, zero, 450, 100, 0.1, 0.40);
  sf_vs_lw_distance_on_w_side = new TH2D("E/p_vs_lw", "E/p_vs_lw", bins, zero, 450, 100, 0.1, 0.40);
  lu_side_distribution = new TH1D("lu_side_distribution", "lu_side_distribution", 50, 0, 400);
  lv_side_distribution = new TH1D("lv_side_distribution", "lv_side_distribution", 50, 0, 450);
  lw_side_distribution = new TH1D("lw_side_distribution", "lw_side_distribution", 50, 0, 450);

  /*
     W_hist_lower = new TH1D("W_lower", "W_lower", bins, zero, w_max);
     Q2_hist_lower = new TH1D("Q2_lower", "Q2_lower", bins, zero, 0.4);
     W_vs_q2_lower = new TH2D("W_vs_q2_lower", "W_vs_q2_lower", bins, zero,
     w_max, bins, zero, 0.4);

     W_hist_upper = new TH1D("W_upper", "W_upper", bins, zero, w_max);
     Q2_hist_upper = new TH1D("Q2_upper", "Q2_upper", bins, 0.4, q2_max);
     W_vs_q2_upper = new TH2D("W_vs_q2_upper", "W_vs_q2_upper", bins, zero,
     w_max, bins, 0.4, q2_max);
   */
  // W_hist_singlePi = new TH1D("W_singlePi", "W_singlePi", bins, zero, w_max);
  // Q2_hist_singlePi = new TH1D("Q2_singlePi", "Q2_singlePi", bins, zero,
  // q2_max);
  // W_vs_q2_singlePi = new TH2D("W_vs_q2_singlePi", "W_vs_q2_singlePi", bins,
  //                             zero, w_max, bins, zero, q2_max);
  /*
     W_hist_lower_singlePi = new TH1D("W_lower_singlePi", "W_lower_singlePi",
     bins, zero, w_max);
     Q2_hist_lower_singlePi = new TH1D("Q2_lower_singlePi", "Q2_lower_singlePi",
     bins, zero, 0.4);
     W_vs_q2_lower_singlePi =
        new TH2D("W_vs_q2_lower_singlePi", "W_vs_q2_lower_singlePi", bins, zero,
     w_max, bins, zero, 0.4);

     W_hist_upper_singlePi = new TH1D("W_upper_singlePi", "W_upper_singlePi",
     bins, zero, w_max);
     Q2_hist_upper_singlePi = new TH1D("Q2_upper_singlePi", "Q2_upper_singlePi",
     bins, 0.4, q2_max);
     W_vs_q2_upper_singlePi =
        new TH2D("W_vs_q2_upper_singlePi", "W_vs_q2_upper_singlePi", bins, zero,
     w_max, bins, 0.4, q2_max);
   */

  makeHists_deltat();
  makeHists_MomVsBeta();
  makeHists_WvsQ2();
  makeHists_MM();
  Make_hist_cc();
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

      hname.append("invariant_mass_2pi");
      htitle.append("invariant_mass_2pi");
      hname.append("_");
      htitle.append(" ");
      hname.append(sec_name[i]);
      htitle.append(sec_name[i]);
      W_hist_2pi[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, p_min, w_max);
      hname.clear();
      htitle.clear();

      hname.append("invariant_mass P#pi+");
      htitle.append("invariant_mass P#pi+");
      hname.append("_");
      htitle.append(" ");
      hname.append(sec_name[i]);
      htitle.append(sec_name[i]);
      W_hist_delta_pp[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, p_min, w_max);
      hname.clear();
      htitle.clear();

      hname.append("invariant_mass P#pi-");
      htitle.append("invariant_mass P#pi-");
      hname.append("_");
      htitle.append(" ");
      hname.append(sec_name[i]);
      htitle.append(sec_name[i]);
      W_hist_delta_zero[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, p_min, w_max);
      hname.clear();
      htitle.clear();

      hname.append("invariant_mass #pi+#pi-");
      htitle.append("invariant_mass #pi+#pi-");
      hname.append("_");
      htitle.append(" ");
      hname.append(sec_name[i]);
      htitle.append(sec_name[i]);
      W_hist_rho[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, p_min, w_max);
      hname.clear();
      htitle.clear();

      hname.append("invariant_mass N#pi+ ");
      htitle.append("invariant_mass N#pi+");
      hname.append("_");
      htitle.append(" ");
      hname.append(sec_name[i]);
      htitle.append(sec_name[i]);
      W_hist_singlepip[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, p_min, w_max);
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
  void Histogram::Fill_WvsmmSQ_ep(double W, double mmSQ, int sec_number) {
    if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
      W_hist_ep[sec_number]->Fill(W);
      W_vs_mmSQ_ep[sec_number][0]->Fill(W, mmSQ);
    }
  }
  void Histogram::Fill_WvsmmSQ_2pi(double W, double W_dpp, double delta_zero_, double rho_, double mmSQ,
                                   int sec_number) {
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

  /*
     if (Q2 <= 0.4) {
      W_vs_q2_lower->Fill(W, Q2);
      W_hist_lower->Fill(W);
      Q2_hist_lower->Fill(Q2);
     } else {
      W_vs_q2_upper->Fill(W, Q2);
      W_hist_upper->Fill(W);
      Q2_hist_upper->Fill(Q2);
     }
   */

  // W and Q^2

  // void Histogram::Fill_WvsQ2_singlePi(double W, double Q2, TLorentzVector
  // *mm)
  // {
  //   W_vs_q2_singlePi->Fill(W, Q2);
  //   W_hist_singlePi->Fill(W);
  //   Q2_hist_singlePi->Fill(Q2);
  //   MM_neutron->Fill(mm->M2());

  /*
     if (Q2 <= 0.4) {
      W_vs_q2_lower_singlePi->Fill(W, Q2);
      W_hist_lower_singlePi->Fill(W);
      Q2_hist_lower_singlePi->Fill(Q2);
     } else {
      W_vs_q2_upper_singlePi->Fill(W, Q2);
      W_hist_upper_singlePi->Fill(W);
      Q2_hist_upper_singlePi->Fill(Q2);
     }
   */
  void Histogram::Make_hist_cc() {
    for (int i = 0; i < cc_num; i++) {
      hname.append("cc_total_");
      htitle.append("cc_total_");
      hname.append(cc_name[i]);
      htitle.append(cc_name[i]);
      cherenkov_total[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 60);
      hname.clear();
      htitle.clear();

      hname.append("cc_ltcc_");
      htitle.append("cc_ltcc_");
      hname.append(cc_name[i]);
      htitle.append(cc_name[i]);
      cherenkov_ltcc[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 20);
      hname.clear();
      htitle.clear();

      hname.append("cc_htcc_");
      htitle.append("cc_htcc_");
      hname.append(cc_name[i]);
      htitle.append(cc_name[i]);
      cherenkov_htcc[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, 0, 60);
      hname.clear();
      htitle.clear();
    }
  }
  void Histogram::Fill_hist_cc_tot(float tot_el) { cherenkov_total[0]->Fill(tot_el); }
  void Histogram::Fill_hist_cc_ltcc(float ltcc_el) { cherenkov_ltcc[0]->Fill(ltcc_el); }
  void Histogram::Fill_hist_cc_htcc(float htcc_el) { cherenkov_htcc[0]->Fill(htcc_el); }
  void Histogram::Fill_hist_cc_tot_pim(float tot_pim) { cherenkov_total[1]->Fill(tot_pim); }
  void Histogram::Fill_hist_cc_ltcc_pim(float ltcc_pim) { cherenkov_ltcc[1]->Fill(ltcc_pim); }
  void Histogram::Fill_hist_cc_htcc_pim(float htcc_pim) { cherenkov_htcc[1]->Fill(htcc_pim); }
  void Histogram::Fill_hist_cc_tot_pip(float tot_pip) { cherenkov_total[2]->Fill(tot_pip); }
  void Histogram::Fill_hist_cc_ltcc_pip(float ltcc_pip) { cherenkov_ltcc[2]->Fill(ltcc_pip); }
  void Histogram::Fill_hist_cc_htcc_pip(float htcc_pip) { cherenkov_htcc[2]->Fill(htcc_pip); }
  //  } else if (cc_num == 1) {

  void Histogram::Write_hist_cc() {
    for (int i = 0; i < cc_num; i++) {
      cherenkov_total[i]->SetXTitle("cc total");
      cherenkov_total[i]->Write();
      delete cherenkov_total[i];

      cherenkov_ltcc[i]->SetXTitle("cc ltcc");
      cherenkov_ltcc[i]->Write();
      delete cherenkov_ltcc[i];

      cherenkov_htcc[i]->SetXTitle("cc htcc");
      cherenkov_htcc[i]->Write();
      delete cherenkov_htcc[i];
    }
  }

  void Histogram::Fill_theta_vs_phi_cm(float th_el, float ph_el) { theta_vs_phi_cm->Fill(th_el, ph_el); }
  void Histogram::Fill_sf_vs_lv(float li, float sf_) { sf_vs_lv_distance_on_v_side->Fill(li, sf_); }
  void Histogram::Fill_sf_vs_lw(float li, float sf_) { sf_vs_lw_distance_on_w_side->Fill(li, sf_); }
  void Histogram::Fill_lu_dist(float li) { lu_side_distribution->Fill(li); }
  void Histogram::Fill_lv_dist(float li) { lv_side_distribution->Fill(li); }
  void Histogram::Fill_lw_dist(float li) { lw_side_distribution->Fill(li); }

  void Histogram::Fill_vertex_vz(float vz) { vertex_vz->Fill(vz); }
  void Histogram::Fill_theta_P(float theta_p, float theta_pip_, float theta_pim_) {
    theta_prot->Fill(theta_p);
    theta_pip->Fill(theta_pip_);
    theta_pim->Fill(theta_pim_);
  }
  void Histogram::Fill_Phi_cm(float Phi_p, float Phi_pip_, float Phi_pim_) {
    Phi_prot->Fill(Phi_p);
    Phi_pip->Fill(Phi_pip_);
    Phi_pim->Fill(Phi_pim_);
  }

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

      W_hist_delta_pp[i]->SetXTitle("W (GeV)");
      W_hist_delta_pp[i]->Write();
      delete W_hist_delta_pp[i];

      W_hist_delta_zero[i]->SetXTitle("W (GeV)");
      W_hist_delta_zero[i]->Write();
      delete W_hist_delta_zero[i];

      W_hist_rho[i]->SetXTitle("W (GeV)");
      W_hist_rho[i]->Write();
      delete W_hist_rho[i];

      W_hist_singlepip[i]->SetXTitle("W (GeV)");
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
    /*
       W_vs_q2_lower->SetXTitle("W (GeV)");
       W_vs_q2_lower->SetYTitle("Q^{2} (GeV^{2})");
       W_vs_q2_lower->SetOption("COLZ");
       W_vs_q2_lower->Write();

       W_hist_lower->SetXTitle("W (GeV)");
       W_hist_lower->Write();

       Q2_hist_lower->SetXTitle("Q^{2} (GeV^{2})");
       Q2_hist_lower->Write();

       W_vs_q2_upper->SetXTitle("W (GeV)");
       W_vs_q2_upper->SetYTitle("Q^{2} (GeV^{2})");
       W_vs_q2_upper->SetOption("COLZ");
       W_vs_q2_upper->Write();

       W_hist_upper->SetXTitle("W (GeV)");
       W_hist_upper->Write();

       Q2_hist_upper->SetXTitle("Q^{2} (GeV^{2})");
       Q2_hist_upper->Write();
     */
    // W_vs_q2_singlePi->SetXTitle("W (GeV)");
    // W_vs_q2_singlePi->SetYTitle("Q^{2} (GeV^{2})");
    // W_vs_q2_singlePi->SetOption("COLZ");
    // W_vs_q2_singlePi->Write();
    //
    // W_hist_singlePi->SetXTitle("W (GeV)");
    // W_hist_singlePi->Write();
    //
    // Q2_hist_singlePi->SetXTitle("Q^{2} (GeV^{2})");
    // Q2_hist_singlePi->Write();
    // MM_neutron->Write();

    theta_prot->SetXTitle("theta_prot_cm");
    theta_prot->Write();
    theta_pip->SetXTitle("theta_pip_cm");
    theta_pip->Write();
    theta_pim->SetXTitle("theta_pim_cm");
    theta_pim->Write();

    Phi_prot->SetXTitle("Phi_prot_cm");
    Phi_prot->Write();
    Phi_pip->SetXTitle("Phi_pip_cm");
    Phi_pip->Write();
    Phi_pim->SetXTitle("Phi_pim_cm");
    Phi_pim->Write();
    /*
       W_vs_q2_lower_singlePi->SetXTitle("W (GeV)");
       W_vs_q2_lower_singlePi->SetYTitle("Q^{2} (GeV^{2})");
       W_vs_q2_lower_singlePi->SetOption("COLZ");
       W_vs_q2_lower_singlePi->Write();

       W_hist_lower_singlePi->SetXTitle("W (GeV)");
       W_hist_lower_singlePi->Write();

       Q2_hist_lower_singlePi->SetXTitle("Q^{2} (GeV^{2})");
       Q2_hist_lower_singlePi->Write();

       W_vs_q2_upper_singlePi->SetXTitle("W (GeV)");
       W_vs_q2_upper_singlePi->SetYTitle("Q^{2} (GeV^{2})");
       W_vs_q2_upper_singlePi->SetOption("COLZ");
       W_vs_q2_upper_singlePi->Write();

       W_hist_upper_singlePi->SetXTitle("W (GeV)");
       W_hist_upper_singlePi->Write();

       Q2_hist_upper_singlePi->SetXTitle("Q^{2} (GeV^{2})");
       Q2_hist_upper_singlePi->Write();
     */
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
          if (m == 1 && e == 1) {
            MM_hist[m][e][i]->Fit("gaus", "", "", -0.05, 0.05);
          }
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
    // for (size_t p = 0; p < particle_num; p++) {
    //   hname.append("delta_t_");
    //   htitle.append("#Deltat ");
    //   hname.append(particle_name[p]);
    //   htitle.append(particle_name[p]);
    //   hname.append("_ctof");
    //   htitle.append(" ctof");
    //   delta_t_hist_ctof[p] = new TH2D(hname.c_str(), htitle.c_str(),
    //   bins,
    //   p_min,
    //                                   p_max, bins, Dt_min, Dt_max);
    //   hname.clear();
    //   htitle.clear();
    // }
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
  }
  // function below is for name and title of histogram vertex
  void Histogram::makeHists_MomVsBeta() {
    for (size_t i = 0; i < with_id_num; i++) {
      hname.append("mom_vs_beta_vertex");
      htitle.append("Momentum vs #beta vertex");
      hname.append("_");
      htitle.append(" ");
      hname.append(id_name[i]);
      htitle.append(id_name[i]);
      momvsbeta_vertex[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, zero, 1.2);
      hname.clear();
      htitle.clear();
    }
    // particle number = 4, 0e 1pi 2P 3 K
    for (size_t p = 0; p < particle_num; p++) {
      for (size_t c = 0; c < charge_num; c++) {
        for (size_t i = 0; i < with_id_num; i++) {
          hname.append("mom_vs_beta_");
          htitle.append("Momentum vs #beta ");
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
          momvsbeta_hist[p][c][i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, p_max, bins, zero, 1.2);
          hname.clear();
          htitle.clear();
        }
      }
    }
  }
  void Histogram::Fill_MomVsBeta_vertex(int pid, int charge, double P, double beta) {
    if (beta != 0) {
      momvsbeta_vertex[0]->Fill(P, beta);
      if (pid == ELECTRON) {
        momvsbeta_vertex[1]->Fill(P, beta);

      } else {
        momvsbeta_vertex[2]->Fill(P, beta);
      }
    }
  }

  void Histogram::Fill_MomVsBeta(int pid, int charge, double P, double beta) {
    int good_ID = 0;
    if (beta != 0) {
      momentum->Fill(P);
      for (size_t p = 0; p < particle_num; p++) {
        switch (p) {
          case 0:
            good_ID = -ELECTRON;
            break;
          case 1:
            good_ID = PIP;
            break;
          case 2:
            good_ID = PROTON;
            break;
          case 3:
            good_ID = KP;
            break;
        }

        /*momvsbeta_hist[p][0][0]->Fill(P, beta);
           if (good_ID == abs(pid)) {
           momvsbeta_hist[p][0][1]->Fill(P, beta);
           } else {
           momvsbeta_hist[p][0][2]->Fill(P, beta);
           }*/

        if (charge == -1) {
          momvsbeta_hist[p][1][0]->Fill(P, beta);
          if (-good_ID == pid) {  // - good_ID thyo paila
            momvsbeta_hist[p][1][1]->Fill(P, beta);
          } else {
            momvsbeta_hist[p][1][2]->Fill(P, beta);
          }
        } else if (charge == 1) {
          momvsbeta_hist[p][0][0]->Fill(P, beta);
          if (good_ID == pid) {
            momvsbeta_hist[p][0][1]->Fill(P, beta);
          } else {
            momvsbeta_hist[p][0][2]->Fill(P, beta);
          }
        }
      }
    }
  }

  void Histogram::Write_MomVsBeta() {
    for (size_t i = 0; i < with_id_num; i++) {
      momvsbeta_vertex[i]->SetXTitle("Momentum (GeV)");
      momvsbeta_vertex[i]->SetYTitle("#beta");
      momvsbeta_vertex[i]->SetOption("COLZ");
      momvsbeta_vertex[i]->Write();
      delete momvsbeta_vertex[i];
    }

    momentum->SetXTitle("Momentum (GeV)");
    momentum->Write();
    for (size_t p = 0; p < particle_num; p++) {
      for (size_t c = 0; c < charge_num; c++) {
        for (size_t i = 0; i < with_id_num; i++) {
          momvsbeta_hist[p][c][i]->SetXTitle("Momentum (GeV)");
          momvsbeta_hist[p][c][i]->SetYTitle("#beta");
          momvsbeta_hist[p][c][i]->SetOption("COLZ");
          momvsbeta_hist[p][c][i]->Write();
          delete momvsbeta_hist[p][c][i];
        }
      }
    }
  }

  void Histogram::Fill_EC(double sf, double momentum) { EC_sampling_fraction->Fill(momentum, sf, int sec); }
  void Histogram::Write_EC() {
    for (int i = 0; i < sec_num; i++) {
      EC_sampling_fraction[i]->SetXTitle("Momentum (GeV)");
      EC_sampling_fraction[i]->SetYTitle("Sampling Fraction");
      EC_sampling_fraction[i]->SetOption("COLZ");
      EC_sampling_fraction[i]->Write();
      delete EC_sampling_fraction[i];
    }
    theta_vs_phi_cm->SetXTitle("theta");
    theta_vs_phi_cm->SetYTitle("phi");
    theta_vs_phi_cm->SetOption("COLZ");
    theta_vs_phi_cm->Write();
    sf_vs_lv_distance_on_v_side->SetXTitle("E/p");
    sf_vs_lv_distance_on_v_side->SetYTitle("distance_on_V_side");
    sf_vs_lv_distance_on_v_side->SetOption("COLZ");
    sf_vs_lv_distance_on_v_side->Write();
    sf_vs_lw_distance_on_w_side->SetXTitle("E/p");
    sf_vs_lw_distance_on_w_side->SetYTitle("distance_on_W_side");
    sf_vs_lw_distance_on_w_side->SetOption("COLZ");
    sf_vs_lw_distance_on_w_side->Write();
    lu_side_distribution->SetXTitle("side_U");
    lu_side_distribution->Write();
    lv_side_distribution->SetXTitle("side_V");
    lv_side_distribution->Write();
    lw_side_distribution->SetXTitle("side_W");
    lw_side_distribution->Write();
    vertex_vz->SetXTitle("vertex_vz");
    vertex_vz->Write();
  }
