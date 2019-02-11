/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
/*  2019-01-31. plot w and mmSq
make different hist for
ctof and ftof
forwardtagger vs forward detector
cut at 2 gev and look */

#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "constants.hpp"
#include "deltat.hpp"
#include "physics.hpp"
#include "reaction.hpp"

class Histogram {
private:
  int bins = 500;
  double p_min = 0.0;
  double p_max = 11.0;
  double Dt_max = 10.0;
  double Dt_min = -Dt_max;
  double q2_max = 8.0;
  double w_max = 5.0;
  double zero = 0.0;

  std::string hname;
  std::string htitle;

  static const short particle_num = 4; // 0-e 1-Pi 2-P 3-K
  std::string particle_name[particle_num] = {"e", "pi", "P", "K"};
  static const short charge_num = 2; // 0-un 1-pos 2-neg
  std::string charge_name[charge_num] = {/*"neutral", */ "positive",
                                         "negative"};
  static const short with_id_num = 3; // 0-without 1-with 2-anti
  std::string id_name[with_id_num] = {"withoutID", "withID", "antiID"};
  static const short sec_num = 6; // 0-without 1-with 2-anti
  std::string sec_name[sec_num] = {"sec_1", "sec_2", "sec_3",
                                   "sec_4", "sec_5", "sec_6"};
  static const short mm_num = 2; // 0 mm 1 mm square
  std::string mm_name[mm_num] = {"mm", "mmSQ"};
  static const short mm_events_num = 8; // 0-ep event 1 2pion ...
  std::string mm_events_name[mm_events_num] = {"elastic",
                                               "2pion",
                                               "Neutron",
                                               "e'p'pi+",
                                               "elastic_wop" /*inclusive*/,
                                               "2pi_wop",
                                               "Neutron_wop",
                                               "e'pi+_wop"};

  // Kinematics

  TH1D *momentum;

  TH2D *W_vs_Q2[sec_num];
  TH1D *W_hist[sec_num];
  TH1D *Q2_hist[sec_num];
  TH1D *MM_hist[mm_num][mm_events_num][sec_num];
  TH2D *W_vs_mmSQ_ep[sec_num];
  TH2D *W_vs_mmSQ_2pi[sec_num];
  TH2D *W_vs_mmSQ_singlepip[sec_num];

  TH1D *W_hist_ep[sec_num];
  TH1D *W_hist_2pi[sec_num];
  TH1D *W_hist_singlepip[sec_num];

  // TH1D *MM_neutron;

  // TH1D *W_hist_lower;
  // TH1D *Q2_hist_lower;
  // TH2D *W_vs_q2_lower;
  //
  // TH1D *W_hist_upper;
  // TH1D *Q2_hist_upper;
  // TH2D *W_vs_q2_upper;
  //
  // TH1D *W_hist_singlePi;
  // TH1D *Q2_hist_singlePi;
  // TH2D *W_vs_q2_singlePi;
  //
  // TH1D *W_hist_lower_singlePi;
  // TH1D *Q2_hist_lower_singlePi;
  // TH2D *W_vs_q2_lower_singlePi;
  //
  // TH1D *W_hist_upper_singlePi;
  // TH1D *Q2_hist_upper_singlePi;
  // TH2D *W_vs_q2_upper_singlePi;

  // EC Sampling Fraction
  TH2D *EC_sampling_fraction;
  // EC Sampling Fraction

  // Mom vs Beta
  TH2D *momvsbeta_hist[particle_num][charge_num][with_id_num];
  TH2D *momvsbeta_vertex[with_id_num];
  // Mom vs Beta

  // Delta T
  TH2D *delta_t_hist[particle_num][charge_num][with_id_num];
  TH2D *delta_t_vertex[with_id_num];
  TH2D *delta_t_ctof_vs_comp;
  TH1D *ctof_comp;
  TH1D *theta_prot;

  // Delta T

public:
  Histogram();
  ~Histogram();
  float mm_lim_max(int mm_number, int mm_events_number);
  float mm_lim_min(int mm_number, int mm_events_number);

  // W and Q^2
  void Fill_ep_mm(double mm, int sec_number);
  void Fill_ep_mmSQ(double mm, int sec_number);
  void Fill_2pion_mm(double mm, int sec_number);
  void Fill_2pion_mmSQ(double mm, int sec_number);
  void Fill_pip_mm(double mm, int sec_number);
  void Fill_pip_mmSQ(double mm, int sec_number);
  void Fill_pim_mm(double mm, int sec_number);
  void Fill_pim_mmSQ(double mm, int sec_number);

  void Fill_MM_wop_e_prime(double mm_1, int sec_number);
  void Fill_MMSQ_wop_e_prime(double mm_1, int sec_number);
  void Fill_MM_wop_2pion(double mm_1, int sec_number);
  void Fill_MMSQ_wop_2pion(double mm_1, int sec_number);
  void Fill_MM_wop_pip(double mm_1, int sec_number);
  void Fill_MMSQ_wop_pip(double mm_1, int sec_number);
  void Fill_MM_wop_pim(double mm_1, int sec_number);
  void Fill_MMSQ_wop_pim(double mm_1, int sec_number);

  // void Fill_WvsQ2_singlePi(double W, double Q2, TLorentzVector *mm);

  void makeHists_WvsQ2();
  void makeHists_MM();
  //  void makeHists_Q2();
  void Fill_WvsmmSQ_ep(double W, double mmSQ, int sec_number);
  void Fill_WvsmmSQ_2pi(double W, double mmSQ, int sec_number);
  void Fill_WvsmmSQ_singlepip(double W, double mmSQ, int sec_number);

  void Fill_WvsQ2(double W, double Q2, int sec_number);
  void Fill_MM_hist(double mm, size_t m, size_t e, int sec_number);
  void Write_WvsQ2();
  void Write_MM_hist();
  void Fill_theta_P(float theta);

  // P and E
  void makeHists_MomVsBeta();
  void Fill_momentum(double P);
  void Fill_MomVsBeta_vertex(int pid, int charge, double P, double beta);
  void Fill_MomVsBeta(int pid, int charge, double P, double beta);
  void Write_MomVsBeta();

  // Delta T
  void makeHists_deltat();
  void Fill_deltat_vertex(int pid, int charge, float dt, float momentum);
  void Fill_deltat_elect(int pid, int charge, float dt, float momentum);
  void Fill_deltat_prot(int pid, int charge, float dt, float momentum);
  void Fill_deltat_pip(int pid, int charge, float dt, float momentum);
  void Fill_deltat_kp(int pid, int charge, float dt, float momentum);

  void Fill_dt_ctof_comp(int ctof_comp, float dt);
  void Fill_ctof_comp(int ctof_comp_1);

  void Write_deltat();

  // EC Sampling Fraction
  void Fill_EC(double etot, double momentum);
  void Write_EC();
};

#endif
