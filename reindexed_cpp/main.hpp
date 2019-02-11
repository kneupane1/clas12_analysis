/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD
#include "TChain.h"
#include "colors.hpp"
#include "constants.hpp"
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
  double tot_energy_ec = 0;
  int sc_d = 0;
  double W = 0;
  double Q2 = 0;
  double sf = 0;
  double per = 0;
  int index = 0;
  int num_pip = 0;
  bool good_e = true;
  int sector;
  Histogram *hist = new Histogram();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (pid->size() == 0 || pid->at(0) != ELECTRON)
      continue;

    per = ((double)current_event / (double)num_of_events);
    if (current_event % 1000 == 0)
      std::cerr << "\t\t" << std::floor(100 * per) << "%\r\r" << std::flush;
    Reaction *event = new Reaction(e_mu);
    event->SetElec(px->at(0), py->at(0), pz->at(0), MASS_E);

    if (event->e_mu_prime().P() != 0)
      hist->Fill_EC(ec_tot_energy->at(0) / event->e_mu_prime().P(),
                    event->e_mu_prime().P());

    Delta_T *dt = new Delta_T(sc_ftof_1b_time->at(0), sc_ftof_1b_path->at(0),
                              sc_ftof_1a_time->at(0), sc_ftof_1a_path->at(0),
                              sc_ftof_2_time->at(0), sc_ftof_2_path->at(0));

    for (int part = 1; part < pid->size(); part++) {
      if (beta->at(part) < 0.02 || p->at(part) < 0.02)
        continue;
      if (dc_sec->at(0) < 7)
        sector = dc_sec->at(0) - 1;

      dt->dt_calc(p->at(part), sc_ftof_1b_time->at(part),
                  sc_ftof_1b_path->at(part), sc_ftof_1a_time->at(part),
                  sc_ftof_1a_path->at(part), sc_ftof_2_time->at(part),
                  sc_ftof_2_path->at(part), sc_ctof_time->at(part),
                  sc_ctof_path->at(part));

      hist->Fill_MomVsBeta_vertex(pid->at(part), charge->at(part), p->at(part),
                                  beta->at(part));

      hist->Fill_MomVsBeta(pid->at(part), charge->at(part), p->at(part),
                           beta->at(part));

      hist->Fill_deltat_vertex(pid->at(0), charge->at(0), dt->dt_E(), p->at(0));

      if (/*event->W() < 2.70 && event->W() > 1.20 &&*/ event->Q2() < 15.0 &&
          event->Q2() > 0.0) {
        if (charge->at(part) == -1) {
          hist->Fill_deltat_elect(pid->at(0), charge->at(0), dt->dt_E(),
                                  p->at(0));
          hist->Fill_deltat_pip(pid->at(part), charge->at(part), dt->dt_Pi(),
                                p->at(part));
          hist->Fill_deltat_kp(pid->at(part), charge->at(part), dt->dt_K(),
                               p->at(part));
          if ((abs(dt->dt_Pi()) < 0.5) /*||
              (dt->dt_Pi() > -4.5 && dt->dt_Pi() < -3.5) /* &&
              (charge->at(part) == -1*/) {
            event->SetPim(px->at(part), py->at(part), pz->at(part), MASS_PIP);
          }
        } else if (charge->at(part) == 1) {
          hist->Fill_deltat_prot(pid->at(part), charge->at(part), dt->dt_P(),
                                 p->at(part));
          hist->Fill_deltat_pip(pid->at(part), charge->at(part), dt->dt_Pi(),
                                p->at(part));
          hist->Fill_deltat_kp(pid->at(part), charge->at(part), dt->dt_K(),
                               p->at(part));
          //    if (charge->at(part) == 1) {
          if ((abs(dt->dt_P()) < 0.5) /* ||
              (dt->dt_P() > -4.5 && dt->dt_P() < -3.7)*/) {
            event->SetProton(px->at(part), py->at(part), pz->at(part), MASS_P);
          } else if ((abs(dt->dt_Pi()) < 0.50) /* ||
                     (dt->dt_Pi() > -4.25 && dt->dt_P() < -3.7 &&
                      event->pip_mu_prime().P() < 1.2)*/) {
            event->SetPip(px->at(part), py->at(part), pz->at(part), MASS_PIP);
          }
        }
        //  }
        if (event->p_mu_prime().Theta() != 0)
          hist->Fill_theta_P(event->p_mu_prime().Theta() * (180 / 3.14));

        dt->dt_calc_1(p->at(part), sc_ctof_time->at(part),
                      sc_ctof_path->at(part));
        if (sc_ctof_component->at(part) > 0) {
          //  std::cout << "dt_ctof" << dt->dt_ctof_P() - dt->dt_P() << "    "
          //  << sc_ctof_component->at(part) << '\n';
          hist->Fill_ctof_comp(sc_ctof_component->at(part));
          hist->Fill_dt_ctof_comp(sc_ctof_component->at(part), dt->dt_ctof_P());
        }
      }
    }
    //  for (int i = 1; i < sector; i++) {
    //  if (event->p_mu_prime().P() != 0) {
    if (/*event->W() < 2.70 && event->W() > 1.20 &&*/ event->Q2() < 15.0 &&
        event->Q2() > 0.0) {
      hist->Fill_WvsQ2(event->W(), event->Q2(), sector);
      //  hist->Fill_WvsmmSQ(event->W(), event->MM2(), sector);
    }
    //}
    //}
    delete dt;

    event->CalcMissMass();

    if (event->elecProtEvent()) {
      hist->Fill_ep_mm(event->MM(), sector);
      hist->Fill_ep_mmSQ(event->MM2(), sector);
      hist->Fill_WvsmmSQ_ep(event->W_ep(), event->MM2(), sector);

    } else if (event->twoPionEvent()) {
      hist->Fill_2pion_mm(event->MM(), sector);
      hist->Fill_2pion_mmSQ(event->MM2(), sector);
      hist->Fill_WvsmmSQ_2pi(event->W_2pi(), event->MM2(), sector);
      //  std::cout << "w_2pi   " << event->W_2pi() << '\n';

      // std::cout << "diff " << event->W_2pi() - event->W() << '\n';
    } else if (event->ProtonPimEvent()) {
      hist->Fill_pip_mm(event->MM(), sector);
      hist->Fill_pip_mmSQ(event->MM2(), sector);
    } else if (event->ProtonPipEvent()) {
      hist->Fill_pim_mm(event->MM(), sector);
      hist->Fill_pim_mmSQ(event->MM2(), sector);
      hist->Fill_WvsmmSQ_singlepip(event->W_singlepip(), event->MM2(), sector);
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
    }
    delete event;
  }

  out->cd();
  hist->Write_EC();
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
