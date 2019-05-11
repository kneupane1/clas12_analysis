/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD
#include <TFile.h>
#include <TLorentzVector.h>
#include <fstream>
#include <vector>
#include "TChain.h"
#include "colors.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "filehandeler.hpp"
#include "histogram.hpp"
#include "physics.hpp"
#include "reaction.hpp"

void datahandeler(std::string fin, std::string fout) {
        double energy = CLAS12_E;
        if (getenv("CLAS12_E") != NULL) energy = atof(getenv("CLAS12_E"));
        TLorentzVector *e_mu = new TLorentzVector(0.0, 0.0, energy, energy);

        TFile *out = new TFile(fout.c_str(), "RECREATE");
        double P;
        bool electron_cuts;
        // Load chain from branch h10
        TChain *chain = filehandeler::addFiles(fin);
        filehandeler::getBranches(chain);

        int num_of_events = (int)chain->GetEntries();
        int total = 0;
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
        int sector;
        float cc_tot;
        float cc_ltcc;
        float cc_htcc;
        float cc_tot_pim;
        float cc_ltcc_pim;
        float cc_htcc_pim;
        float cc_tot_pip;
        float cc_ltcc_pip;
        float cc_htcc_pip;
        int n_pip = 0;
        int n_pim = 0;

        Histogram *hist = new Histogram();

        for (int current_event = 0; current_event < num_of_events; current_event++) {
                chain->GetEntry(current_event);
                if (pid->size() == 0 || /*pid->at(0) != ELECTRON*/ charge->at(0) >= 0) continue;

                per = ((double)current_event / (double)num_of_events);
                if (current_event % 1000 == 0) std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;
                Reaction *event = new Reaction(e_mu);
                event->SetElec(px->at(0), py->at(0), pz->at(0), MASS_E);

                Cuts *e_cuts = new Cuts();
                good_e = e_cuts->electron_cuts(status->at(0), charge->at(0), event->e_mu_prime().P(),
                                               (ec_tot_energy->at(0) / event->e_mu_prime().P()), vz->at(0), chi2pid->at(0));
                if (good_e == false) continue;
                if (event->e_mu_prime().P() != 0)
                        hist->Fill_EC(ec_tot_energy->at(0) / event->e_mu_prime().P(), event->e_mu_prime().P());

                Delta_T *dt = new Delta_T(sc_ftof_1b_time->at(0), sc_ftof_1b_path->at(0), sc_ftof_1a_time->at(0),
                                          sc_ftof_1a_path->at(0), sc_ftof_2_time->at(0), sc_ftof_2_path->at(0));
                if (pid->at(0) == ELECTRON) {
                        cc_tot = cc_nphe_tot->at(0);
                        if (cc_tot >= 0) {
                                hist->Fill_hist_cc_tot(cc_tot);
                        }
                        cc_ltcc = cc_ltcc_nphe->at(0);
                        if (cc_ltcc >= 0) {
                                hist->Fill_hist_cc_ltcc(cc_ltcc);
                        }
                        cc_htcc = cc_htcc_nphe->at(0);
                        if (cc_htcc >= 0) {
                                hist->Fill_hist_cc_htcc(cc_htcc);
                        }
                }
                n_pip = 0;
                n_pim = 0;
                for (int part = 1; part < pid->size(); part++) {
                        // if (pid->at(part) == PIP) {
                        //   n_pip = n_pip + 1;
                        // }
                        // if (n_pip > 1) continue;
                        // if (pid->at(part) == PIM) {
                        //   n_pim = n_pim + 1;
                        // }
                        // if (n_pim > 1) continue;

                        if (beta->at(part) < 0.02 || p->at(part) < 0.02) continue;
                        if (dc_sec->at(0) < 7) {
                                sector = dc_sec->at(0) - 1;
                        }

                        dt->dt_calc(p->at(part), sc_ftof_1b_time->at(part), sc_ftof_1b_path->at(part), sc_ftof_1a_time->at(part),
                                    sc_ftof_1a_path->at(part), sc_ftof_2_time->at(part), sc_ftof_2_path->at(part), sc_ctof_time->at(part),
                                    sc_ctof_path->at(part));

                        hist->Fill_MomVsBeta_vertex(pid->at(part), charge->at(part), p->at(part), beta->at(part));

                        hist->Fill_MomVsBeta(pid->at(part), charge->at(part), p->at(part), beta->at(part));

                        hist->Fill_deltat_vertex(pid->at(0), charge->at(0), dt->dt_E(), p->at(0));

                        if (/*event->W() < 1.40 &&  event->W() > 1.20 &&*/ event->Q2() < 15.0 && event->Q2() > 0.0) {
                                if (charge->at(part) == -1) {
                                        hist->Fill_deltat_elect(pid->at(0), charge->at(0), dt->dt_E(), p->at(0));
                                        hist->Fill_deltat_pip(pid->at(part), charge->at(part), dt->dt_Pi(), p->at(part));
                                        hist->Fill_deltat_kp(pid->at(part), charge->at(part), dt->dt_K(), p->at(part));
                                        if ((abs(dt->dt_Pi()) < 0.5) ||
                                            (dt->dt_Pi() > -4.5 && dt->dt_Pi() < -3.5) /* &&
                                                                                                                                                                                    (charge->at(part) == -1*/) {
<<<<<<< HEAD
            event->SetPim(px->at(part), py->at(part), pz->at(part), MASS_PIP);
            good_pim = e_cuts->pim_cuts(status->at(part), charge->at(part), event->pim_mu_prime().P(), pid->at(part),
                                        chi2pid->at(part));
            if (pid->at(part) == PIM) {
              cc_tot_pim = cc_nphe_tot->at(part);
              if (cc_tot_pim >= 0) {
                hist->Fill_hist_cc_tot_pim(cc_tot_pim);
              }
              cc_ltcc_pim = cc_ltcc_nphe->at(part);
              if (cc_ltcc_pim >= 0) {
                hist->Fill_hist_cc_ltcc_pim(cc_ltcc_pim);
              }
              cc_htcc_pim = cc_htcc_nphe->at(part);
              if (cc_htcc_pim >= 0) {
                hist->Fill_hist_cc_htcc_pim(cc_htcc_pim);
              }
            }
          }
        } else if (charge->at(part) == 1) {
          hist->Fill_deltat_prot(pid->at(part), charge->at(part), dt->dt_P(), p->at(part));
          hist->Fill_deltat_pip(pid->at(part), charge->at(part), dt->dt_Pi(), p->at(part));
          hist->Fill_deltat_kp(pid->at(part), charge->at(part), dt->dt_K(), p->at(part));
          //    if (charge->at(part) == 1) {
          if ((abs(dt->dt_P()) < 0.5) || (dt->dt_P() > -4.5 && dt->dt_P() < -3.7)) {
            event->SetProton(px->at(part), py->at(part), pz->at(part), MASS_P);
            good_p = e_cuts->proton_cuts(status->at(part), charge->at(part), event->p_mu_prime().P(), pid->at(part),
                                         chi2pid->at(part));
          } else if ((dt->dt_Pi() > -0.10) && (dt->dt_Pi() < 0.50) ||
                     ((dt->dt_Pi() > -4.25 && dt->dt_P() < -3.7) && event->pip_mu_prime().P() < 0.7)) {
            event->SetPip(px->at(part), py->at(part), pz->at(part), MASS_PIP);
            good_pip = e_cuts->pip_cuts(status->at(part), charge->at(part), event->pip_mu_prime().P(), pid->at(part),
                                        chi2pid->at(part));
            // if (pid->at(part) == 2212 && charge->at(part) > 0) {
            if (pid->at(part) == PIP) {
              cc_tot_pip = cc_nphe_tot->at(part);
              if (cc_tot_pip >= 0) {
                hist->Fill_hist_cc_tot_pip(cc_tot_pip);
              }
              cc_ltcc_pip = cc_ltcc_nphe->at(part);
              if (cc_ltcc_pip >= 0) {
                hist->Fill_hist_cc_ltcc_pip(cc_ltcc_pip);
              }
              cc_htcc_pip = cc_htcc_nphe->at(part);
              if (cc_htcc_pip >= 0) {
                hist->Fill_hist_cc_htcc_pip(cc_htcc_pip);
              }
            }
          }
        } else
          event->SetOther(px->at(part), py->at(part), pz->at(part), MASS_N, pid->at(part));
      }
    }
    // if (good_pim == false) continue;
    // if (good_p == false) continue;
    // if (good_pim == false) continue;

    // dt->dt_calc_1(p->at(part), sc_ctof_time->at(part), sc_ctof_path->at(part));
    // if (sc_ctof_component->at(part) > 0) {
    //         //  std::cout << "dt_ctof" << dt->dt_ctof_P() - dt->dt_P() << "    "
    //  << sc_ctof_component->at(part) << '\n';

    hist->Fill_sf_vs_lu((ec_ecin_lu->at(0) + ec_ecout_lu->at(0) + ec_pcal_lu->at(0)) / 3,
                        (ec_tot_energy->at(0) / event->e_mu_prime().P()));
    hist->Fill_sf_vs_lv((ec_ecin_lv->at(0) + ec_ecout_lv->at(0) + ec_pcal_lv->at(0)) / 3,
                        (ec_tot_energy->at(0) / event->e_mu_prime().P()));
    hist->Fill_sf_vs_lw((ec_ecin_lw->at(0) + ec_ecout_lw->at(0) + ec_pcal_lw->at(0)) / 3,
                        (ec_tot_energy->at(0) / event->e_mu_prime().P()));
    hist->Fill_lu_dist((ec_ecin_lu->at(0) + ec_ecout_lu->at(0) + ec_pcal_lu->at(0)) / 3);
    hist->Fill_lv_dist((ec_ecin_lv->at(0) + ec_ecout_lv->at(0) + ec_pcal_lv->at(0)) / 3);
    hist->Fill_lw_dist((ec_ecin_lw->at(0) + ec_ecout_lw->at(0) + ec_pcal_lw->at(0)) / 3);

    hist->Fill_vertex_vz(vz->at(0));

    // if (event->twoPionEvent()) {

    //  for (int i = 1; i < sector; i++) {
    // if (event->p_mu_prime().P() != 0) {
    hist->Fill_WvsQ2(event->W(), event->Q2(), sector);
    //}
    //}
    delete dt;
    if (/*event->W() < 1.40 &&  event->W() > 1.20 &&*/ event->Q2() < 15.0 && event->Q2() > 0.0) {
      event->CalcMissMass();
      event->AlphaCalc();
      if (event->twoPionEvent()) {
        if (event->p_mu_prime_cm().Theta() > 0) {
          hist->Fill_theta_P(event->p_mu_prime_cm().Theta() * (180 / PI), event->pip_mu_prime_cm().Theta() * (180 / PI),
                             event->pim_mu_prime_cm().Theta() * (180 / PI));
=======
                                                event->SetPim(px->at(part), py->at(part), pz->at(part), MASS_PIP);
                                                good_pim = e_cuts->pim_cuts(status->at(part), charge->at(part), event->pim_mu_prime().P(), pid->at(part),
                                                                            chi2pid->at(part));
                                                if (pid->at(part) == PIM) {
                                                        cc_tot_pim = cc_nphe_tot->at(part);
                                                        if (cc_tot_pim >= 0) {
                                                                hist->Fill_hist_cc_tot_pim(cc_tot_pim);
                                                        }
                                                        cc_ltcc_pim = cc_ltcc_nphe->at(part);
                                                        if (cc_ltcc_pim >= 0) {
                                                                hist->Fill_hist_cc_ltcc_pim(cc_ltcc_pim);
                                                        }
                                                        cc_htcc_pim = cc_htcc_nphe->at(part);
                                                        if (cc_htcc_pim >= 0) {
                                                                hist->Fill_hist_cc_htcc_pim(cc_htcc_pim);
                                                        }
                                                }
                                        }
                                } else if (charge->at(part) == 1) {
                                        hist->Fill_deltat_prot(pid->at(part), charge->at(part), dt->dt_P(), p->at(part));
                                        hist->Fill_deltat_pip(pid->at(part), charge->at(part), dt->dt_Pi(), p->at(part));
                                        hist->Fill_deltat_kp(pid->at(part), charge->at(part), dt->dt_K(), p->at(part));
                                        //    if (charge->at(part) == 1) {
                                        if ((abs(dt->dt_P()) < 0.5) || (dt->dt_P() > -4.5 && dt->dt_P() < -3.7)) {
                                                event->SetProton(px->at(part), py->at(part), pz->at(part), MASS_P);
                                                good_p = e_cuts->proton_cuts(status->at(part), charge->at(part), event->p_mu_prime().P(), pid->at(part),
                                                                             chi2pid->at(part));
                                        } else if ((abs(dt->dt_Pi()) < 0.50) ||
                                                   (dt->dt_Pi() > -4.25 && dt->dt_P() < -3.7) /*&&
                                                                                                              event->pip_mu_prime().P() < 1.2)*/) {
                                                event->SetPip(px->at(part), py->at(part), pz->at(part), MASS_PIP);
                                                good_pip = e_cuts->pip_cuts(status->at(part), charge->at(part), event->pip_mu_prime().P(), pid->at(part),
                                                                            chi2pid->at(part));
                                                // if (pid->at(part) == 2212 && charge->at(part) > 0) {
                                                if (pid->at(part) == PIP) {
                                                        cc_tot_pip = cc_nphe_tot->at(part);
                                                        if (cc_tot_pip >= 0) {
                                                                hist->Fill_hist_cc_tot_pip(cc_tot_pip);
                                                        }
                                                        cc_ltcc_pip = cc_ltcc_nphe->at(part);
                                                        if (cc_ltcc_pip >= 0) {
                                                                hist->Fill_hist_cc_ltcc_pip(cc_ltcc_pip);
                                                        }
                                                        cc_htcc_pip = cc_htcc_nphe->at(part);
                                                        if (cc_htcc_pip >= 0) {
                                                                hist->Fill_hist_cc_htcc_pip(cc_htcc_pip);
                                                        }
                                                }
                                        }
                                } else
                                        event->SetOther(px->at(part), py->at(part), pz->at(part), MASS_N, pid->at(part));
                        }
                }
                if (good_pip == false) continue;
                if (good_p == false) continue;
                if (good_pim == false) continue;

                // dt->dt_calc_1(p->at(part), sc_ctof_time->at(part), sc_ctof_path->at(part));
                // if (sc_ctof_component->at(part) > 0) {
                //         //  std::cout << "dt_ctof" << dt->dt_ctof_P() - dt->dt_P() << "    "
                //  << sc_ctof_component->at(part) << '\n';

                hist->Fill_sf_vs_lu((ec_ecin_lu->at(0) + ec_ecout_lu->at(0) + ec_pcal_lu->at(0)) / 3,
                                    (ec_tot_energy->at(0) / event->e_mu_prime().P()));
                hist->Fill_sf_vs_lv((ec_ecin_lv->at(0) + ec_ecout_lv->at(0) + ec_pcal_lv->at(0)) / 3,
                                    (ec_tot_energy->at(0) / event->e_mu_prime().P()));
                hist->Fill_sf_vs_lw((ec_ecin_lw->at(0) + ec_ecout_lw->at(0) + ec_pcal_lw->at(0)) / 3,
                                    (ec_tot_energy->at(0) / event->e_mu_prime().P()));
                hist->Fill_lu_dist((ec_ecin_lu->at(0) + ec_ecout_lu->at(0) + ec_pcal_lu->at(0)) / 3);
                hist->Fill_lv_dist((ec_ecin_lv->at(0) + ec_ecout_lv->at(0) + ec_pcal_lv->at(0)) / 3);
                hist->Fill_lw_dist((ec_ecin_lw->at(0) + ec_ecout_lw->at(0) + ec_pcal_lw->at(0)) / 3);

                hist->Fill_vertex_vz(vz->at(0));

                // if (event->twoPionEvent()) {

                //  for (int i = 1; i < sector; i++) {
                // if (event->p_mu_prime().P() != 0) {
                hist->Fill_WvsQ2(event->W(), event->Q2(), sector);
                //}
                //}
                delete dt;
                if (/*event->W() < 1.40 &&  event->W() > 1.20 &&*/ event->Q2() < 15.0 && event->Q2() > 0.0) {
                        event->CalcMissMass();
                        event->AlphaCalc();
                        if (event->twoPionEvent()) {
                                if (event->p_mu_prime_cm().Theta() > 0) {
                                        hist->Fill_theta_P(event->p_mu_prime_cm().Theta() * (180 / PI), event->pip_mu_prime_cm().Theta() * (180 / PI),
                                                           event->pim_mu_prime_cm().Theta() * (180 / PI));
                                }

                                hist->Fill_Phi_cm(event->p_mu_prime_cm().Phi() * (180 / PI), event->pip_mu_prime_cm().Phi() * (180 / PI),
                                                  event->pim_mu_prime_cm().Phi() * (180 / PI));

                                //{
                                // std::cout << "electron " << event->q_cm().M() << '\n';

                                // n += 1;
                                // std::cout
                                //<< "  pip_pim_alpha_cm  " << event->alpha_ppip_pipim()
                                //<< "  pip_theta_cm  "
                                //<< event->pip_mu_prime_cm().Theta() * (180 / PI)
                                //<< "   pim_theta_cm  "
                                //  << event->pim_mu_prime_cm().Theta() * (180 / PI) << "  " << n
                                //  << '\n';
                        }
                        if (event->elecProtEvent()) {
                                hist->Fill_ep_mm(event->MM(), sector);
                                hist->Fill_ep_mmSQ(event->MM2(), sector);
                                if (event->MM2() < 0.1 && event->MM2() > -0.1) {
                                        hist->Fill_WvsmmSQ_ep(event->W_ep(), event->MM2(), sector);
                                } else {
                                        hist->Fill_WvsmmSQ_anti_ep(event->W_ep(), event->MM2(), sector);
                                }

                        } else if (event->twoPionEvent()) {
                                hist->Fill_2pion_mm(event->MM(), sector);
                                hist->Fill_2pion_mmSQ(event->MM2(), sector);
                                if (event->MM2() < 0.1 && event->MM2() > -0.1) {
                                        hist->Fill_WvsmmSQ_2pi(event->W_2pi(), event->W_delta_pp(), event->W_delta_zero(), event->W_rho(),
                                                               event->MM2(), sector);
                                } else {
                                        hist->Fill_WvsmmSQ_anti_2pi(event->W_2pi(), event->W_delta_pp(), event->W_delta_zero(), event->W_rho(),
                                                                    event->MM2(), sector);
                                }
                                //  std::cout << "w_2pi   " << event->W_2pi() << '\n';

                                // std::cout << "diff " << event->W_2pi() - event->W() << '\n';
                        } else if (event->ProtonPimEvent()) {
                                hist->Fill_pip_mm(event->MM(), sector);
                                hist->Fill_pip_mmSQ(event->MM2(), sector);

                        } else if (event->ProtonPipEvent()) {
                                hist->Fill_pim_mm(event->MM(), sector);
                                hist->Fill_pim_mmSQ(event->MM2(), sector);
                        }

                        event->CalcMissMass_wop();

                        if (event->elecWopEvent()) {
                                hist->Fill_MM_wop_e_prime(event->MM_wop(), sector);
                                hist->Fill_MMSQ_wop_e_prime(event->MM2_wop(), sector);
                        } else if (event->twoPionWopEvent()) {
                                hist->Fill_MM_wop_2pion(event->MM_wop(), sector);
                                hist->Fill_MMSQ_wop_2pion(event->MM2_wop(), sector);

                        } else if (event->WopPimEvent()) {
                                hist->Fill_MM_wop_pip(event->MM_wop(), sector);
                                hist->Fill_MMSQ_wop_pip(event->MM2_wop(), sector);
                        } else if (event->WopPipEvent()) {
                                hist->Fill_MM_wop_pim(event->MM_wop(), sector);
                                hist->Fill_MMSQ_wop_pim(event->MM2_wop(), sector);
                                if (event->MM2_wop() < 1.05 && event->MM2_wop() > 0.7) {
                                        hist->Fill_WvsmmSQ_singlepip(event->W_singlepip(), event->MM2_wop(), sector);
                                } else {
                                        hist->Fill_WvsmmSQ_anti_singlepip(event->W_singlepip(), event->MM2_wop(), sector);
                                }
                        }
                }
                delete e_cuts;
                delete event;
>>>>>>> b7677dc757f2ac927e69665971bc730557a76473
        }

        out->cd();
        TDirectory *CC_EC = out->mkdir("ccEc");
        CC_EC->cd();
        hist->Write_EC();
        hist->Write_hist_cc();

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
