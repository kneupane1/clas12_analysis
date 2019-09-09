/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD
#include "TChain.h"
#include "colors.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "filehandeler.hpp"
#include "histogram.hpp"
#include "physics.hpp"
#include "reaction.hpp"
#include <TFile.h>
#include <TLorentzVector.h>
#include <fstream>
#include <vector>

void datahandeler(std::string fin, std::string fout) {
  double energy = CLAS12_E;
  if (getenv("CLAS12_E") != NULL)
    energy = atof(getenv("CLAS12_E"));
  TLorentzVector *e_mu = new TLorentzVector(0.0, 0.0, energy, energy);

  TFile *out = new TFile(fout.c_str(), "RECREATE");
  double P;
  bool electron_cuts;
  // Load chain from branch h10
  TChain *chain = filehandeler::addFiles(fin);
  filehandeler::getBranches(chain);

  int num_of_events = (int)chain->GetEntries();
  int total = 0;

  int pid_of_part = 0;
  float mom_at_part = 0.0;
  int charge_at_part = 0.0;

  float e_pim_mom_diff = 0;
  float e_pim_mom_diff_x = 0;
  float e_pim_mom_diff_y = 0;
  float e_pim_mom_diff_z = 0;
  float p_pip_mom_diff = 0;
  float p_pip_mom_diff_x = 0;
  float p_pip_mom_diff_y = 0;
  float p_pip_mom_diff_z = 0;
  double tot_energy_ec = 0;

  int sc_d = 0;
  double W = 0;
  double Q2 = 0;
  double sf = 0;
  double per = 0;
  int index = 0;
  int num_pip = 0;
  bool good_e = false;
  bool good_p = false;
  bool good_pip = false;
  bool good_pim = false;
  bool good_hadron_ctof_pim = false;
  bool good_hadron_ctof_pip = false;
  bool good_hadron_ctof_P = false;

  int sector;
  // float cc_tot;
  // float cc_ltcc;
  // float cc_htcc;
  // float cc_tot_pim;
  // float cc_ltcc_pim;
  // float cc_htcc_pim;
  // float cc_tot_pip;
  // float cc_ltcc_pip;
  // float cc_htcc_pip;

  Histogram *hist = new Histogram();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);

    if (pid->size() == 0 || pid->at(0) != ELECTRON /* charge->at(0) >= 0*/)
      continue; // cut # 1

    per = ((double)current_event / (double)num_of_events);
    if (current_event % 1000 == 0)
      std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;
    Reaction *event = new Reaction(e_mu);
    event->SetElec(px->at(0), py->at(0), pz->at(0), MASS_E);
    if (dc_sec->at(0) < 7) {
      sector = dc_sec->at(0) - 1;
    }

    hist->Fill_hist_PCAL_without_FID_CUT(ec_pcal_x->at(0), ec_pcal_y->at(0));
    hist->Fill_hist_DC_without_FID_CUT(dc_r1_x->at(0), dc_r1_y->at(0),
                                       dc_r2_x->at(0), dc_r2_y->at(0),
                                       dc_r3_x->at(0), dc_r3_y->at(0));

    Cuts *e_cuts = new Cuts();
    good_e = e_cuts->electron_cuts(
        status->at(0), charge->at(0),
        (ec_tot_energy->at(0) / event->e_mu_prime().P()), vz->at(0),
        chi2pid->at(0), event->e_mu_prime().P(),
        event->e_mu_prime().Theta() * 180 / PI,
        event->e_mu_prime().Phi() * 180 / PI, sector, ec_pcal_x->at(0),
        ec_pcal_y->at(0), dc_r1_x->at(0), dc_r1_y->at(0), dc_r2_x->at(0),
        dc_r2_y->at(0), dc_r3_x->at(0), dc_r3_y->at(0));

    if (good_e) { // cut # 2

      if (event->e_mu_prime().P() != 0) {
        hist->Fill_EC_sampling_fraction(
            event->e_mu_prime().P(),
            (ec_tot_energy->at(0) / event->e_mu_prime().P()), sector,
            mc_weight);
        hist->Fill_PCAL_VS_ECAL(
            ec_pcal_energy->at(0),
            (ec_ecin_energy->at(0) + ec_ecout_energy->at(0)), sector,
            mc_weight);
        hist->Fill_hist_PCAL_FID_CUT(ec_pcal_x->at(0), ec_pcal_y->at(0));
        hist->Fill_hist_DC_FID_CUT(dc_r1_x->at(0), dc_r1_y->at(0),
                                   dc_r2_x->at(0), dc_r2_y->at(0),
                                   dc_r3_x->at(0), dc_r3_y->at(0));
      }

      Delta_T *dt = new Delta_T(sc_ftof_1b_time->at(0), sc_ftof_1b_path->at(0),
                                sc_ftof_1a_time->at(0), sc_ftof_1a_path->at(0),
                                sc_ftof_2_time->at(0), sc_ftof_2_path->at(0));
      // if (pid->at(0) == ELECTRON) {
      //   cc_tot = cc_nphe_tot->at(0);
      //   if (cc_tot >= 0) {
      //     hist->Fill_hist_cc_tot(cc_tot);
      //   }
      //   cc_ltcc = cc_ltcc_nphe->at(0);
      //   if (cc_ltcc >= 0) {
      //     hist->Fill_hist_cc_ltcc(cc_ltcc);
      //   }
      //   cc_htcc = cc_htcc_nphe->at(0);
      //   if (cc_htcc >= 0) {
      //     hist->Fill_hist_cc_htcc(cc_htcc);
      //   }
      // }

      for (int part = 1; part < pid->size(); part++) {
        pid_of_part = pid->at(part);
        mom_at_part = p->at(part);
        charge_at_part = charge->at(part);

        hist->Fill_pid_size_fn(pid->size(), pid->at(part));

        if (beta->at(part) < 0.02 || p->at(part) < 0.02)
          continue; // cut # 3

        dt->dt_calc(p->at(part), sc_ftof_1b_time->at(part),
                    sc_ftof_1b_path->at(part), sc_ftof_1a_time->at(part),
                    sc_ftof_1a_path->at(part), sc_ftof_2_time->at(part),
                    sc_ftof_2_path->at(part), sc_ctof_time->at(part),
                    sc_ctof_path->at(part));

        dt->dt_calc_1(p->at(part), sc_ctof_time->at(part),
                      sc_ctof_path->at(part));

        hist->Fill_MomVsBeta_vertex(pid->at(part), charge->at(part),
                                    p->at(part), beta->at(part));

        hist->Fill_MomVsBeta(pid->at(part), charge->at(part), p->at(part),
                             beta->at(part));

        hist->Fill_deltat_vertex(pid->at(0), charge->at(0), dt->dt_E(),
                                 p->at(0));

        if (/*event->W() < 1.40 &&  event->W() > 1.20 &&*/ event->Q2() < 15.0 &&
            event->Q2() > 0.0) { // cut # 5

          if (charge->at(part) == -1) { // cut # 6

            if (abs(dt->dt_Pi()) < 10.1) { // cut # 7
              hist->Fill_deltat_pion(pid->at(part), charge->at(part),
                                     dt->dt_Pi(), p->at(part));
            }
            if (abs(dt->dt_ctof_Pi()) < 10.1) {
              hist->Fill_ctof_pim_without_cut_hist(
                  pid->at(part), charge->at(part), dt->dt_ctof_Pi(),
                  p->at(part));

              if (abs(dt->dt_ctof_Pi()) < 0.5) {

                hist->Fill_ctof_pim_with_cut_hist(
                    pid->at(part), charge->at(part), dt->dt_ctof_Pi(),
                    p->at(part));
              }
            }
            if ((abs(dt->dt_Pi()) < 0.5 || abs(dt->dt_ctof_Pi()) < 0.5) &&
                (pid->at(part) == PIM)) { // cut # 8

              event->SetPim(px->at(part), py->at(part), pz->at(part), MASS_PIP);
              good_pim = e_cuts->pim_cuts(status->at(part), charge->at(part),
                                          event->pim_mu_prime().P(),
                                          pid->at(part), chi2pid->at(part));
              good_hadron_ctof_pim = e_cuts->hadron_cuts_ctof(
                  status->at(part), charge->at(part), event->pim_mu_prime().P(),
                  pid->at(part), chi2pid->at(part));

              //  }
            }

            //  if (abs(dt->dt_E()) < 10.1) {  // cut # 7
            hist->Fill_deltat_elect(pid->at(0), charge->at(0), dt->dt_E(),
                                    p->at(0));
            //}

            // if (abs(dt->dt_ctof_E()) < 10.1) {
            // hist->Fill_ctof_e_without_cut_hist(pid->at(0), charge->at(0),
            //                                 dt->dt_ctof_E(), p->at(0));
            // if (abs(dt->dt_ctof_E()) < 0.5)
            // hist->Fill_ctof_e_with_cut_hist(pid->at(0), charge->at(0),
            //                              dt->dt_ctof_E(), p->at(0));
            //}
            // if (abs(dt->dt_K()) < 10.1) {  //// cut # 7
            hist->Fill_deltat_kp(pid->at(part), charge->at(part), dt->dt_K(),
                                 p->at(part));
            //}
            // if (abs(dt->dt_ctof_K()) < 10.1) {
            hist->Fill_ctof_km_without_cut_hist(pid->at(part), charge->at(part),
                                                dt->dt_ctof_K(), p->at(part));
            if (abs(dt->dt_ctof_K()) < 0.5) {
              hist->Fill_ctof_km_with_cut_hist(pid->at(part), charge->at(part),
                                               dt->dt_ctof_K(), p->at(part));
            }

            // } else if (((abs(dt->dt_K()) < 0.5) ||
            //             ((abs(dt->dt_ctof_K()) < 0.5) || (dt->dt_ctof_K() >
            //             -4.5 && dt->dt_ctof_K() < -3.5))) &&
            //            (pid->at(part) == KM)) {          // cut # 8
            //         event->SetKm(px->at(part), py->at(part), pz->at(part),
            //         MASS_KM); if (pid->at(part) == PIM) {
            //                 cc_tot_pim = cc_nphe_tot->at(part);
            //                 if (cc_tot_pim >= 0) {
            //                         hist->Fill_hist_cc_tot_pim(cc_tot_pim);
            //                 }
            //                 cc_ltcc_pim = cc_ltcc_nphe->at(part);
            //                 if (cc_ltcc_pim >= 0) {
            //                         hist->Fill_hist_cc_ltcc_pim(cc_ltcc_pim);
            //                 }
            //                 cc_htcc_pim = cc_htcc_nphe->at(part);
            //                 if (cc_htcc_pim >= 0) {
            //                         hist->Fill_hist_cc_htcc_pim(cc_htcc_pim);
            //                }
            //         }
          } else if (charge->at(part) == 1) { // cut # 6
            if (abs(dt->dt_P()) < 10.0) {     // cut # 7
              hist->Fill_deltat_prot(pid->at(part), charge->at(part),
                                     dt->dt_P(), p->at(part));
            }
            if (abs(dt->dt_ctof_P()) < 10.0) {
              hist->Fill_ctof_P_without_cut_hist(pid->at(part),
                                                 charge->at(part),
                                                 dt->dt_ctof_P(), p->at(part));
              if (abs(dt->dt_ctof_P()) < 0.5) {

                hist->Fill_ctof_P_with_cut_hist(pid->at(part), charge->at(part),
                                                dt->dt_ctof_P(), p->at(part));
              }
            }
            if (abs(dt->dt_Pi()) < 10.0) { // cut # 7
              hist->Fill_deltat_pion(pid->at(part), charge->at(part),
                                     dt->dt_Pi(), p->at(part));
            }
            if (abs(dt->dt_ctof_Pi()) < 10.0) {
              hist->Fill_ctof_pip_without_cut_hist(
                  pid->at(part), charge->at(part), dt->dt_ctof_Pi(),
                  p->at(part));
              if ((dt->dt_ctof_Pi() < 0.5) && (dt->dt_ctof_Pi() > -0.30)) {
                //(dt->dt_ctof_Pi() >= (0.67 * p->at(part) - 4.5)))
                hist->Fill_ctof_pip_with_cut_hist(
                    pid->at(part), charge->at(part), dt->dt_ctof_Pi(),
                    p->at(part));
              }
            }
            if ((abs(dt->dt_P()) < 0.5 || abs(dt->dt_ctof_P()) < 0.5)
                //(dt->dt_ctof_P() > (1.92* p->at(part) * p->at(part)) - 3.673 *
                // p->at(part) +2)))
                && (pid->at(part) == PROTON)) { // cut 9

              event->SetProton(px->at(part), py->at(part), pz->at(part),
                               MASS_P);
              good_p = e_cuts->proton_cuts(status->at(part), charge->at(part),
                                           event->p_mu_prime().P(),
                                           pid->at(part), chi2pid->at(part));
              good_hadron_ctof_P = e_cuts->hadron_cuts_ctof(
                  status->at(part), charge->at(part), event->p_mu_prime().P(),
                  pid->at(part), chi2pid->at(part));

            } else if ((abs(dt->dt_Pi()) < 0.50 ||
                        (dt->dt_ctof_Pi() < 0.5 && dt->dt_ctof_Pi() > -0.3)) &&
                       pid->at(part) == PIP) { // cut 9
              event->SetPip(px->at(part), py->at(part), pz->at(part), MASS_PIP);
              good_pip = e_cuts->pip_cuts(status->at(part), charge->at(part),
                                          event->pip_mu_prime().P(),
                                          pid->at(part), chi2pid->at(part));
              good_hadron_ctof_pip = e_cuts->hadron_cuts_ctof(
                  status->at(part), charge->at(part), event->pip_mu_prime().P(),
                  pid->at(part), chi2pid->at(part));
            }

            // if (abs(dt->dt_K()) < 10.1) {  // cut # 7
            hist->Fill_deltat_kp(pid->at(part), charge->at(part), dt->dt_K(),
                                 p->at(part));
            //}
            // if (abs(dt->dt_ctof_K()) < 10.1) {
            hist->Fill_ctof_kp_without_cut_hist(pid->at(part), charge->at(part),
                                                dt->dt_ctof_K(), p->at(part));
            if ((((dt->dt_K() > -0.50) && (dt->dt_K() < 0.50)) ||
                                                     (dt->dt_ctof_K() < 0.5 && dt->dt_ctof_K() > -0.5)) /*&&
                                                                                                            (pid->at(part) == KP)*/) {
              hist->Fill_ctof_kp_with_cut_hist(pid->at(part), charge->at(part),
                                               dt->dt_ctof_K(), p->at(part));
            }
            //}

            // if (pid->at(part) == 2212 && charge->at(part) > 0) {
            // if (pid->at(part) == PIP) {
            //   cc_tot_pip = cc_nphe_tot->at(part);
            //   if (cc_tot_pip >= 0) {
            //     hist->Fill_hist_cc_tot_pip(cc_tot_pip);
            //   }
            //   cc_ltcc_pip = cc_ltcc_nphe->at(part);
            //   if (cc_ltcc_pip >= 0) {
            //     hist->Fill_hist_cc_ltcc_pip(cc_ltcc_pip);
            //   }
            //   cc_htcc_pip = cc_htcc_nphe->at(part);
            //   if (cc_htcc_pip >= 0) {
            //     hist->Fill_hist_cc_htcc_pip(cc_htcc_pip);
            //}
            //}
          } else // cut # 6
            event->SetOther(px->at(part), py->at(part), pz->at(part), MASS_N,
                            pid->at(part));
        } // cut # 5
      }

      // hist->Fill_lu_dist((ec_ecin_lu->at(0) + ec_ecout_lu->at(0) +
      // ec_pcal_lu->at(0)) / 3); hist->Fill_lv_dist((ec_ecin_lv->at(0) +
      // ec_ecout_lv->at(0) + ec_pcal_lv->at(0)) / 3);
      // hist->Fill_lw_dist((ec_ecin_lw->at(0) + ec_ecout_lw->at(0) +
      // ec_pcal_lw->at(0)) / 3);

      hist->Fill_vertex_vz(vz->at(0));
      // hist->Fill_pid_size_fn(pid->size(), (pid_check_el));
      // if (event->twoPionEvent()) {

      //  for (int i = 1; i < sector; i++) {
      // if (event->p_mu_prime().P() != 0) {
      hist->Fill_WvsQ2(event->W(), event->Q2(), sector, mc_weight);

      hist->Fill_W_vs_Q2_all_sec(event->W(), event->Q2(), mc_weight);
      if ((good_pip || good_hadron_ctof_pip) &&
          event->pip_mu_prime().Theta() != 0 &&
          event->pip_mu_prime().P() != 0) {

        hist->Fill_hist_mass_vs_q2_pip(
            event->W(), event->pip_mu_prime().Mag(),
            (event->pip_mu_prime().Theta() * (180 / PI)), event->Q2(),
            mc_weight);

        hist->Fill_mom_pip(event->pip_mu_prime().P(),
                           event->pip_mu_prime().Pz());
      }
      if ((good_pim || good_hadron_ctof_pim) &&
          event->pim_mu_prime().Theta() != 0 &&
          event->pim_mu_prime().P() != 0) {
        hist->Fill_hist_mass_vs_q2_pim(
            event->W(), event->pim_mu_prime().Mag(),
            (event->pim_mu_prime().Theta() * (180 / PI)), event->Q2(),
            mc_weight);

        hist->Fill_mom_pim(event->pim_mu_prime().P(),
                           event->pim_mu_prime().Pz());
      }

      if ((good_p || good_hadron_ctof_P) && event->p_mu_prime().Theta() != 0 &&
          event->p_mu_prime().P() != 0) {

        hist->Fill_hist_mass_vs_q2_prot(
            event->W(), event->p_mu_prime().Mag(),
            (event->p_mu_prime().Theta() * (180 / PI)), event->Q2(), mc_weight);
        hist->Fill_mom_p(event->p_mu_prime().P(), event->p_mu_prime().Pz());
      }

      if (event->elecPimEvent() && good_pim) {

        e_pim_mom_diff = event->e_mu_prime().P() - event->pim_mu_prime().P();
        e_pim_mom_diff_x =
            event->e_mu_prime().Px() - event->pim_mu_prime().Px();
        e_pim_mom_diff_y =
            event->e_mu_prime().Py() - event->pim_mu_prime().Py();
        e_pim_mom_diff_z =
            event->e_mu_prime().Pz() - event->pim_mu_prime().Pz();
        hist->Fill_mom_diff_e_pim(e_pim_mom_diff, e_pim_mom_diff_x,
                                  e_pim_mom_diff_y, e_pim_mom_diff_z);

        //}
      }
      if (event->ProtonPipEvent() && (good_p || good_hadron_ctof_P) &&
          (good_pip || good_hadron_ctof_pip)) {
        p_pip_mom_diff = event->p_mu_prime().P() - event->pip_mu_prime().P();
        p_pip_mom_diff_x =
            event->p_mu_prime().Px() - event->pip_mu_prime().Px();
        p_pip_mom_diff_y =
            event->p_mu_prime().Py() - event->pip_mu_prime().Py();
        p_pip_mom_diff_z =
            event->p_mu_prime().Pz() - event->pip_mu_prime().Pz();
        hist->Fill_mom_diff_p_pip(p_pip_mom_diff, p_pip_mom_diff_x,
                                  p_pip_mom_diff_y, p_pip_mom_diff_z);
      }
      if (/*event->W() < 1.40 &&  event->W() > 1.20 &&*/ event->Q2() < 15.0 &&
          event->Q2() > 0.0 && (good_p || good_hadron_ctof_P) &&
          (good_pip || good_hadron_ctof_pip) &&
          (good_pim || good_hadron_ctof_pim)) { // cut # 10 11 12
        event->CalcMissMass();
        // event->AlphaCalc();      // for now i am not using alpha
        // function from reaction class
        if (event->elecProtEvent()) {
          //  if (good_p || good_hadron_ctof_P) {
          hist->Fill_ep_mm(event->MM(), sector, mc_weight);
          hist->Fill_ep_mmSQ(event->MM2(), sector, mc_weight);
          if (event->MM2() < 0.2 && event->MM2() > -0.1) {
            hist->Fill_WvsmmSQ_ep(event->W_ep(), event->MM2(), sector,
                                  mc_weight);

          } else {
            hist->Fill_WvsmmSQ_anti_ep(event->W_ep(), event->MM2(), sector,
                                       mc_weight);
          }
          //}
        } else if (event->twoPionEvent() /* && (good_p || good_hadron_ctof_P) &&
                   (good_pip == true || good_hadron_ctof_pip == true) &&
                   (good_pim == true || good_hadron_ctof_pim == true)*/) { //
          //     cut # 13

          hist->Fill_2pion_mm(event->MM(), sector, mc_weight);
          hist->Fill_2pion_mmSQ(event->MM2(), sector, mc_weight);
          // if (0.18< ec_tot_energy->at(0) / event->e_mu_prime().P() &&
          // ec_tot_energy->at(0) / event->e_mu_prime().P()< 0.28) {
          //
          //         hist->Fill_2pion_mmSQ(event->MM2(), sector,
          //         mc_weight);
          // }
          if (event->MM2() < 0.03 && event->MM2() > -0.03) { // cut # 14

            // if (abs(dt->dt_P()) < 0.6)
            hist->Fill_deltat_prot_after(pid_of_part, charge_at_part,
                                         dt->dt_P(), dt->dt_Pi(), mom_at_part);
            // if (dt->dt_ctof_P() < 0.5 && dt->dt_ctof_P() > -0.5)
            hist->Fill_deltat_ctof_prot_after(pid_of_part, charge_at_part,
                                              dt->dt_ctof_P(), dt->dt_ctof_Pi(),
                                              mom_at_part);

            hist->Fill_WvsmmSQ_2pi(event->W_2pi(), event->W_delta_pp(),
                                   event->W_delta_zero(), event->W_rho(),
                                   event->MM2(), sector, mc_weight);
            hist->Fill_W_2pi_all_sec(event->W_2pi(), event->W_delta_pp(),
                                     event->W_delta_zero(), event->W_rho(),
                                     mc_weight);

            hist->Fill_theta_P_inv_mass(
                event->W_rho(), (event->p_mu_prime_cm().Theta() * (180 / PI)),
                mc_weight);

            hist->Fill_theta_pim_inv_mass(
                event->W_delta_pp(),
                (event->pim_mu_prime_cm().Theta() * (180 / PI)), mc_weight);
            hist->Fill_theta_pip_inv_mass(
                event->W_delta_zero(),
                (event->pip_mu_prime_cm().Theta() * (180 / PI)), mc_weight);

            hist->Fill_theta_P_lab_inv_mass(
                event->W_rho(), (event->p_mu_prime().Theta() * (180 / PI)),
                mc_weight);
            hist->Fill_theta_pim_lab_inv_mass(
                event->W_delta_pp(),
                (event->pim_mu_prime().Theta() * (180 / PI)), mc_weight);
            hist->Fill_theta_pip_lab_inv_mass(
                event->W_delta_zero(),
                (event->pip_mu_prime().Theta() * (180 / PI)), mc_weight);
          } else {
            hist->Fill_WvsmmSQ_anti_2pi(event->W_2pi(), event->W_delta_pp(),
                                        event->W_delta_zero(), event->W_rho(),
                                        event->MM2(), sector, mc_weight);
          }
        } else if (event->ProtonPimEvent()) {
          hist->Fill_pip_mm(event->MM(), sector, mc_weight);
          hist->Fill_pip_mmSQ(event->MM2(), sector, mc_weight);

        } else if (event->ProtonPipEvent()) {
          hist->Fill_pim_mm(event->MM(), sector, mc_weight);
          hist->Fill_pim_mmSQ(event->MM2(), sector, mc_weight);
        }

        event->CalcMissMass_wop();

        if (event->elecWopEvent()) {
          hist->Fill_MM_wop_e_prime(event->MM_wop(), sector, mc_weight);
          hist->Fill_MMSQ_wop_e_prime(event->MM2_wop(), sector, mc_weight);
        } else if (event->twoPionWopEvent()) {
          hist->Fill_MM_wop_2pion(event->MM_wop(), sector, mc_weight);
          hist->Fill_MMSQ_wop_2pion(event->MM2_wop(), sector, mc_weight);

        } else if (event->WopPimEvent()) {
          hist->Fill_MM_wop_pip(event->MM_wop(), sector, mc_weight);
          hist->Fill_MMSQ_wop_pip(event->MM2_wop(), sector, mc_weight);
        } else if (event->WopPipEvent()) {
          hist->Fill_MM_wop_pim(event->MM_wop(), sector, mc_weight);
          hist->Fill_MMSQ_wop_pim(event->MM2_wop(), sector, mc_weight);
          if (event->MM_wop() < 1.05 && event->MM2_wop() > 0.65) {
            hist->Fill_deltat_after_mmsq_cut_Xpip(pid_of_part, charge_at_part,
                                                  dt->dt_P(), dt->dt_Pi(),
                                                  mom_at_part);
            hist->Fill_deltat_ctof_after_mmsq_cut_Xpip(
                pid_of_part, charge_at_part, dt->dt_ctof_P(), dt->dt_ctof_Pi(),
                mom_at_part);

            hist->Fill_WvsmmSQ_singlepip(event->W_singlepip(), event->MM2_wop(),
                                         sector, mc_weight);
            hist->Fill_W_hist_Xpip_all_sec(event->W_singlepip(), mc_weight);

          } else {
            hist->Fill_WvsmmSQ_anti_singlepip(
                event->W_singlepip(), event->MM2_wop(), sector, mc_weight);
          }
        }
      }

      //}
      //}
      delete dt;

    } // cut # 2

    delete e_cuts;
    delete event;
  }

  out->cd();
  TDirectory *CC_EC = out->mkdir("ccEc");
  CC_EC->cd();
  hist->Write_EC();
  // hist->Write_hist_cc();

  TDirectory *MM = out->mkdir("missingMass");
  MM->cd();
  hist->Write_MM_hist();

  TDirectory *wvsq2 = out->mkdir("wvsq2");
  wvsq2->cd();
  hist->Write_WvsQ2();

  TDirectory *mom_vs_beta = out->mkdir("mom_vs_beta");
  mom_vs_beta->cd();
  hist->Write_MomVsBeta();

  TDirectory *deltat_ftof = out->mkdir("deltat_ftof");
  deltat_ftof->cd();
  hist->Write_deltat();

  out->Close();
  chain->Reset();
  std::cerr << "\nErrors: " << total << "\t" << std::endl;
  delete hist;
}
#endif
