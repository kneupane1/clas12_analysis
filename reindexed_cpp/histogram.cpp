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
  theta_prot =
      new TH1D("theta_distribution", "theta_dist-prot", bins, zero, 90);

  ctof_comp = new TH1D("ctof_component", "ctof_component", bins, zero, 55);

  delta_t_ctof_vs_comp =
      new TH2D("delta_t_ctof_vs_comp", "delta_t_ctof_vs_comp", 50, zero, 50,
               bins, Dt_min, Dt_max);
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
  EC_sampling_fraction =
      new TH2D("EC_sampling_fraction", "EC_sampling_fraction", bins, p_min,
               p_max, bins, zero, 1.0);

  makeHists_deltat();
  makeHists_MomVsBeta();
  makeHists_WvsQ2();
  makeHists_MM();
}

Histogram::~Histogram() {}
// W and Q^2

float Histogram::mm_lim_min(int mm_number, int mm_events_number) {

  if (mm_number == 0 && mm_events_number < 4) {
    return -4.0;
  } else if (mm_number == 0 && mm_events_number > 4) {
    return -4.0;
  } else if (mm_number == 0 && mm_events_number == 4) {
    return 0;
  } else if (mm_number == 1 && mm_events_number < 4) {
    return -4;
  } else if (mm_number == 1 && mm_events_number > 4) {
    return -4;
  } else if (mm_number == 1 && mm_events_number == 4) {
    return -1.5;
  } else {
    return -5;
  }
}
float Histogram::mm_lim_max(int mm_number, int mm_events_number) {

  if (mm_number == 0 && mm_events_number < 4) {
    return 4.0;
  } else if (mm_number == 0 && mm_events_number > 4) {
    return 5.0;
  } else if (mm_number == 0 && mm_events_number == 4) {
    return 5.9;
  } else if (mm_number == 1 && mm_events_number < 4) {
    return 6;
  } else if (mm_number == 1 && mm_events_number > 4) {
    return 16.0;
  } else if (mm_number == 1 && mm_events_number == 4) {
    return 20.5;
  } else {
    return 10;
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
    W_vs_Q2[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min, w_max,
                          bins, p_min, q2_max);
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

    hname.append("Q2_hist");
    htitle.append("Q2_hist");
    hname.append("_");
    htitle.append(" ");
    hname.append(sec_name[i]);
    htitle.append(sec_name[i]);
    Q2_hist[i] = new TH1D(hname.c_str(), htitle.c_str(), bins, p_min, q2_max);
    hname.clear();
    htitle.clear();
  }
}

void Histogram::Fill_WvsQ2(double W, double Q2, int sec_number) {

  if (sec_number == sec_number && sec_number >= 0 && sec_number < 7) {
    W_vs_Q2[sec_number]->Fill(W, Q2);
    W_hist[sec_number]->Fill(W);
    Q2_hist[sec_number]->Fill(Q2);
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

// void Histogram::Fill_WvsQ2_singlePi(double W, double Q2, TLorentzVector *mm)
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

void Histogram::Fill_dt_ctof_comp(int ctof_comp, float dt) {
  delta_t_ctof_vs_comp->Fill(ctof_comp, dt);
}
void Histogram::Fill_ctof_comp(int ctof_comp_1) {
  ctof_comp->Fill(ctof_comp_1);
}
void Histogram::Fill_theta_P(float theta) { theta_prot->Fill(theta); }

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

    Q2_hist[i]->SetXTitle("Q^{2} (GeV^{2})");
    Q2_hist[i]->Write();
    delete Q2_hist[i];
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

  delta_t_ctof_vs_comp->SetXTitle("ctof_component");
  delta_t_ctof_vs_comp->SetYTitle("dt_ctof");
  delta_t_ctof_vs_comp->SetOption("COLZ");
  delta_t_ctof_vs_comp->Write();

  ctof_comp->SetXTitle("ctof_comp");
  ctof_comp->Write();
  theta_prot->SetXTitle("theta_prot");
  theta_prot->Write();

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
            new TH1D(hname.c_str(), htitle.c_str(), bins,
                     Histogram::mm_lim_min(m, e), Histogram::mm_lim_max(m, e));
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
    delta_t_vertex[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min,
                                 p_max, bins, Dt_min, Dt_max);
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
        delta_t_hist[p][c][i] = new TH2D(hname.c_str(), htitle.c_str(), bins,
                                         p_min, p_max, bins, Dt_min, Dt_max);
        hname.clear();
        htitle.clear();
      }
    }
  }
}

void Histogram::Fill_deltat_vertex(int pid, int charge, float dt,
                                   float momentum) {
  delta_t_vertex[0]->Fill(momentum, dt);
  if (pid == ELECTRON) {
    delta_t_vertex[1]->Fill(momentum, dt);
  } else {
    delta_t_vertex[2]->Fill(momentum, dt);
  }
}
void Histogram::Fill_deltat_elect(int pid, int charge, float dt,
                                  float momentum) {
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

void Histogram::Fill_deltat_prot(int pid, int charge, float dt,
                                 float momentum) {
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
    momvsbeta_vertex[i] = new TH2D(hname.c_str(), htitle.c_str(), bins, p_min,
                                   p_max, bins, zero, 1.2);
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
        momvsbeta_hist[p][c][i] = new TH2D(hname.c_str(), htitle.c_str(), bins,
                                           p_min, p_max, bins, zero, 1.2);
        hname.clear();
        htitle.clear();
      }
    }
  }
}
void Histogram::Fill_MomVsBeta_vertex(int pid, int charge, double P,
                                      double beta) {
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
        if (-good_ID == pid) { // - good_ID thyo paila
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

void Histogram::Fill_EC(double sf, double momentum) {
  EC_sampling_fraction->Fill(momentum, sf);
}
void Histogram::Write_EC() {
  EC_sampling_fraction->SetXTitle("Momentum (GeV)");
  EC_sampling_fraction->SetYTitle("Sampling Fraction");
  EC_sampling_fraction->SetOption("COLZ");
  EC_sampling_fraction->Write();
}
