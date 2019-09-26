/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

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
float p_min = 0.0;
float p_max = 11.0;
float Dt_max = 10.0;
float Dt_min = -Dt_max;
float q2_max = 8.0;
float w_max = 5.0;
float zero = 0.0;

std::string hname;
std::string htitle;

static const short particle_num = 4;    // 0-e 1-Pi 2-P 3-K
std::string particle_name[particle_num] = {"e", "pi", "P", "K"};
static const short charge_num = 2;    // 0-un 1-pos 2-neg
std::string charge_name[charge_num] = { /*"neutral", */ "positive", "negative"};
static const short with_id_num = 3;    // 0-without 1-with 2-anti
std::string id_name[with_id_num] = {"withoutID", "withID", "antiID"};
static const short sec_num = 6;    // 0-without 1-with 2-anti
std::string sec_name[sec_num] = {"sec_1", "sec_2", "sec_3", "sec_4", "sec_5", "sec_6"};
static const short mm_num = 2;    // 0 mm 1 mm square
std::string mm_name[mm_num] = {"mm", "mmSQ"};
static const short mm_events_num = 8;    // 0-ep event 1 2pion ...
std::string mm_events_name[mm_events_num] = {
        "e(p,p'X)e'",     "e(p,p'pi+pi-X)e'", "e(p,p'pi-X)e'", "e(p,p'pi+X)e'", "e(p,X)e'" /*inclusive*/,
        "e(p,pi+pi-X)e'", "e(p,pi-X)e'",      "e(p,pi+X)e'"
};

static const short cc_num = 3;    // 0-without 1-with 2-anti
std::string cc_name[cc_num] = {"ele", " pi-", "pi+" /*"total", "htcc", "ltcc"*/};
static const short cut_y_n = 2;
std::string cut_name[cut_y_n] = {"with_cut", "with_anti_cut"};
static const short cut_without_cut_num = 2;
std::string cut_without_cut_name[cut_y_n] = {"with_cut", "with_out_cut"};
static const short slice_num = 20;
std::string slice_name[slice_num] = {"1-1.3",   "1.3-1.6", "1.6-1.9", "1.9-2.2", "2.2-2.5", "2.5-2.8", "2.8-3.1",
                                     "3.1-3.4", "3.4-3.7", "3.7-4.0", "4.0-4.3", "4.3-4.6", "4.6-4.9", "4.9-5.2",
                                     "5.2-5.5", "5.5-5.8", "5.8-6.1", "6.1-6.4", "6.4-6.7", "6.7-7.0"};

static const short w_range_num = 4;
std::string w_range_name[w_range_num] = {" W<2.0 ", " 2.0 <W< 2.5 ", " 2.5 <W< 3.0 ", " 3.0<Q2<3.5 "};

static const short q2_range_num = 3;
std::string q2_range_name[q2_range_num] = {" Q2<2.0 ", " 2.0<Q2<5.0 ", " 5.0<Q2<10.0 "};
// Kinematics

TH1D *pid_size1;
TH1D *pid_check;
TH1D *mom_diff_e_pim;
TH1D *mom_diff_p_pip;
TH1D *mom_diff_e_pim_x;
TH1D *mom_diff_p_pip_x;
TH1D *mom_diff_e_pim_y;
TH1D *mom_diff_p_pip_y;
TH1D *mom_diff_e_pim_z;
TH1D *mom_diff_p_pip_z;

TH1D *mom_e;
TH1D *mom_p;
TH1D *mom_pip;
TH1D *mom_pim;
TH1D *mom_e_z;
TH1D *mom_p_z;
TH1D *mom_pip_z;
TH1D *mom_pim_z;

TH1D *momentum;
TH1D *cherenkov_total[cc_num];
TH1D *cherenkov_ltcc[cc_num];
TH1D *cherenkov_htcc[cc_num];

TH2D *Deltat_P_after_mmsq_cut;
TH2D *Deltat_pip_after_mmsq_cut;

TH2D *Deltat_ctof_P_after_mmsq_cut;
TH2D *Deltat_ctof_pip_after_mmsq_cut;

TH2D *Deltat_P_after_mmsq_cut_singlepip;
TH2D *Deltat_pip_after_mmsq_cut_singlepip;

TH2D *Deltat_ctof_P_after_mmsq_cut_singlepip;
TH2D *Deltat_ctof_pip_after_mmsq_cut_singlepip;

TH1D *W_thrown;
TH2D *W_vs_Q2_thrown;
TH2D *W_vs_Q2_all_sec;

TH2D *Prot_mass_w_vs_Q2[w_range_num][q2_range_num];
TH2D *Prot_theta_lab_vs_W[w_range_num][q2_range_num];
TH2D *Pip_mass_w_vs_Q2[w_range_num][q2_range_num];
TH2D *Pip_theta_lab_vs_W[w_range_num][q2_range_num];
TH2D *Pim_theta_lab_vs_W[w_range_num][q2_range_num];

TH2D *W_vs_Q2[sec_num];
TH1D *W_hist[sec_num];
TH1D *Q2_hist[sec_num];
TH1D *MM_hist[mm_num][mm_events_num][sec_num];
TH2D *W_vs_mmSQ_ep[sec_num][cut_y_n];
TH2D *W_vs_mmSQ_2pi[sec_num][cut_y_n];
TH2D *W_vs_mmSQ_singlepip[sec_num][cut_y_n];

TH1D *W_hist_ep[sec_num];
TH1D *W_hist_2pi[sec_num];
TH1D *W_hist_twopi_all_sec[w_range_num][q2_range_num];
TH1D *inv_mass_P_pip_all_sec[w_range_num][q2_range_num];
TH1D *inv_mass_P_pim_all_sec[w_range_num][q2_range_num];
TH1D *inv_mass_pip_pim_all_sec[w_range_num][q2_range_num];
TH1D *W_hist_delta_pp[sec_num];
TH1D *W_hist_delta_zero[sec_num];
TH1D *W_hist_rho[sec_num];

TH2D *theta_P_vs_mass_pip_pim[w_range_num][q2_range_num];
TH2D *theta_pim_vs_mass_Ppip[w_range_num][q2_range_num];
TH2D *theta_pip_vs_mass_Ppim[w_range_num][q2_range_num];
TH2D *theta_P_lab_vs_mass_pip_pim[w_range_num][q2_range_num];
TH2D *theta_pim_lab_vs_mass_Ppip[w_range_num][q2_range_num];
TH2D *theta_pip_lab_vs_mass_Ppim[w_range_num][q2_range_num];

TH2D *theta_pim_vs_mass_Ppip_thrown[w_range_num][q2_range_num];
TH2D *theta_pip_vs_mass_Ppim_thrown[w_range_num][q2_range_num];
TH2D *theta_P_vs_mass_pip_pim_thrown[w_range_num][q2_range_num];

TH2D *theta_pim_lab_vs_mass_Ppip_thrown[w_range_num][q2_range_num];
TH2D *theta_pip_lab_vs_mass_Ppim_thrown[w_range_num][q2_range_num];
TH2D *theta_P_lab_vs_mass_pip_pim_thrown[w_range_num][q2_range_num];
TH1D *W_hist_singlepip[sec_num];
TH1D *W_hist_Xpip_all_sec;

TH1D *W_hist_twopi_thrown[w_range_num][q2_range_num];
TH1D *inv_mass_P_pip_thrown[w_range_num][q2_range_num];
TH1D *inv_mass_P_pim_thrown[w_range_num][q2_range_num];
TH1D *inv_mass_pip_pim_thrown[w_range_num][q2_range_num];

// angles in cm frame
TH1D *theta_P_cm[w_range_num][q2_range_num];
TH1D *theta_pip_cm[w_range_num][q2_range_num];
TH1D *theta_pim_cm[w_range_num][q2_range_num];
TH1D *phi_P_cm[w_range_num][q2_range_num];
TH1D *phi_pip_cm[w_range_num][q2_range_num];
TH1D *phi_pim_cm[w_range_num][q2_range_num];
TH1D *alpha_P_cm[w_range_num][q2_range_num];
TH1D *alpha_pip_cm[w_range_num][q2_range_num];
TH1D *alpha_pim_cm[w_range_num][q2_range_num];

TH1D *theta_P_cm_thrown[w_range_num][q2_range_num];
TH1D *theta_pip_cm_thrown[w_range_num][q2_range_num];
TH1D *theta_pim_cm_thrown[w_range_num][q2_range_num];
TH1D *phi_P_cm_thrown[w_range_num][q2_range_num];
TH1D *phi_pip_cm_thrown[w_range_num][q2_range_num];
TH1D *phi_pim_cm_thrown[w_range_num][q2_range_num];
TH1D *alpha_P_cm_thrown[w_range_num][q2_range_num];
TH1D *alpha_pip_cm_thrown[w_range_num][q2_range_num];
TH1D *alpha_pim_cm_thrown[w_range_num][q2_range_num];
// EC Sampling Fraction
TH1D *sf_[sec_num][slice_num];
TH2D *EC_sampling_fraction[sec_num];
TH2D *PCAL_VS_ECAL[sec_num];
TH2D *PCAL_FID_CUT[cut_y_n];
TH2D *DCr1_FID_CUT[cut_y_n];
TH2D *DCr2_FID_CUT[cut_y_n];
TH2D *DCr3_FID_CUT[cut_y_n];

// Mom vs Beta
TH2D *momvsbeta_hist[particle_num][charge_num][with_id_num];
TH2D *momvsbeta_vertex[with_id_num];
// Mom vs Beta

// Delta T
TH2D *delta_t_hist[particle_num][charge_num][with_id_num];
TH2D *delta_t_vertex[with_id_num];
TH2D *delta_t_hist_ctof[particle_num][charge_num][cut_without_cut_num];

TH1D *lu_side_distribution;
TH1D *lv_side_distribution;
TH1D *lw_side_distribution;

TH1D *vertex_vz;

// Delta T

public:
Histogram();
~Histogram();
float mm_lim_max(int mm_number, int mm_events_number);
float mm_lim_min(int mm_number, int mm_events_number);

float w_hist_lim_max(int w_num);

// W and Q^2
void Fill_ep_mm(float mm, int sec_number, float weight);
void Fill_ep_mmSQ(float mm, int sec_number, float weight);
void Fill_2pion_mm(float mm, int sec_number, float weight);
void Fill_2pion_mmSQ(float mm, int sec_number, float weight);
void Fill_pip_mm(float mm, int sec_number, float weight);
void Fill_pip_mmSQ(float mm, int sec_number, float weight);
void Fill_pim_mm(float mm, int sec_number, float weight);
void Fill_pim_mmSQ(float mm, int sec_number, float weight);

void Fill_MM_wop_e_prime(float mm_1, int sec_number, float weight);
void Fill_MMSQ_wop_e_prime(float mm_1, int sec_number, float weight);
void Fill_MM_wop_2pion(float mm_1, int sec_number, float weight);
void Fill_MMSQ_wop_2pion(float mm_1, int sec_number, float weight);
void Fill_MM_wop_pip(float mm_1, int sec_number, float weight);
void Fill_MMSQ_wop_pip(float mm_1, int sec_number, float weight);
void Fill_MM_wop_pim(float mm_1, int sec_number, float weight);
void Fill_MMSQ_wop_pim(float mm_1, int sec_number, float weight);

void Fill_W_2pi_all_sec(float W, float q2, float W_2pi, float W_dpp, float delta_zero_, float rho_, float weight);
void Fill_W_hist_Xpip_all_sec(float W, float weight);

void Fill_W_2pi_thrown(float W, float q2, float W_2pi, float W_dpp, float delta_zero_, float rho_, float weight);

void makeHists_WvsQ2();
void makeHists_Angles();
void makeHists_MM();
void makeHists_EC_sf();
void makeHists_sf_slices();
void makeHists_pcal_fid_cuts();

void Fill_W_vs_Q2_all_sec(float w, float q2, float wt);

void Fill_Angles(float W, float q2, float theta_p, float theta_pip, float theta_pim,
                 float phi_p, float phi_pip, float phi_pim,
                 float alpha_p, float alpha_pip, float alpha_pim,float weight);
void Fill_Angles_thrown(float W, float q2, float theta_p, float theta_pip, float theta_pim,
                        float phi_p, float phi_pip, float phi_pim,
                        float alpha_p, float alpha_pip, float alpha_pim,float weight);

// void Fill_WvsQ2_range(float W, float Q2, float weight);

void Fill_hist_mass_vs_q2_prot(float w, float m_p, float theta, float q2, float wt);
void Fill_hist_mass_vs_q2_pip(float w, float m_pip, float theta, float q2, float wt);
void Fill_hist_mass_vs_q2_pim(float w, float m_pim, float theta, float q2, float wt);

//  void makeHists_Q2();
void Fill_WvsmmSQ_ep(float W, float mmSQ, int sec_number, float weight);
void Fill_WvsmmSQ_2pi(float W, float W_dpp, float delta_zero_, float rho_, float mmSQ, int sec_number, float weight);
void Fill_WvsmmSQ_singlepip(float W, float mmSQ, int sec_number, float weight);
void Fill_WvsmmSQ_anti_ep(float W, float mmSQ, int sec_number, float weight);
void Fill_WvsmmSQ_anti_2pi(float W, float W_dpp, float delta_zero_, float rho_, float mmSQ, int sec_number,
                           float weight);
void Fill_WvsmmSQ_anti_singlepip(float W, float mmSQ, int sec_number, float weight);
void Fill_WvsQ2(float W, float Q2, int sec_number, float weight);
void Fill_W_vs_Q2_thrown(float w, float q2, float wt);
void Fill_MM_hist(float mm, size_t m, size_t e, int sec_number, float weight);

void Write_WvsQ2();
void Write_Angles();
void Write_MM_hist();
void Fill_theta_P(float theta_p, float theta_pip_, float theta_pim_);
void Fill_Phi_cm(float Phi_p, float Phi_pip_, float Phi_pim_);
// P and E
void makeHists_MomVsBeta();
void Fill_momentum(float P);

void Make_hist_cc();
void Fill_hist_cc_tot(float tot_el /*, float ltcc_el, float htcc_el*/);
void Fill_hist_cc_ltcc(/*float tot_el, */ float ltcc_el /*, float htcc_el*/);
void Fill_hist_cc_htcc(/*float tot_el, float ltcc_el, */ float htcc_el);
void Fill_hist_cc_tot_pim(float tot_pim);
void Fill_hist_cc_ltcc_pim(float ltcc_pim /*, float htcc_el*/);
void Fill_hist_cc_htcc_pim(float htcc_pim);
void Fill_hist_cc_tot_pip(float tot_pip);
void Fill_hist_cc_ltcc_pip(float ltcc_pip /*, float htcc_el*/);
void Fill_hist_cc_htcc_pip(float htcc_pip);
void Write_hist_cc();

void Fill_MomVsBeta_vertex(int pid, int charge, float P, float beta);
void Fill_MomVsBeta(int pid, int charge, float P, float beta);
void Write_MomVsBeta();

// Delta T
void makeHists_deltat();
void Fill_deltat_vertex(int pid, int charge, float dt, float momentum);
void Fill_deltat_elect(int pid, int charge, float dt, float momentum);
void Fill_deltat_prot(int pid, int charge, float dt, float momentum);
void Fill_deltat_pion(int pid, int charge, float dt, float momentum);
void Fill_deltat_kp(int pid, int charge, float dt, float momentum);

void Fill_ctof_e_without_cut_hist(int pid, int charge, float dt_ctof, float momentum);
void Fill_ctof_P_without_cut_hist(int pid, int charge, float dt_ctof, float momentum);
void Fill_ctof_pip_without_cut_hist(int pid, int charge, float dt_ctof, float momentum);
void Fill_ctof_pim_without_cut_hist(int pid, int charge, float dt_ctof, float momentum);
void Fill_ctof_kp_without_cut_hist(int pid, int charge, float dt_ctof, float momentum);
void Fill_ctof_km_without_cut_hist(int pid, int charge, float dt_ctof, float momentum);

void Fill_ctof_e_with_cut_hist(int pid, int charge, float dt_ctof, float momentum);
void Fill_ctof_P_with_cut_hist(int pid, int charge, float dt_ctof, float momentum);
void Fill_ctof_pip_with_cut_hist(int pid, int charge, float dt_ctof, float momentum);
void Fill_ctof_pim_with_cut_hist(int pid, int charge, float dt_ctof, float momentum);
void Fill_ctof_kp_with_cut_hist(int pid, int charge, float dt_ctof, float momentum);
void Fill_ctof_km_with_cut_hist(int pid, int charge, float dt_ctof, float momentum);

void Fill_deltat_prot_after(int pid, int charge, float dt_p, float dt_pip, float momentum);
void Fill_deltat_ctof_prot_after(int pid, int charge, float dt_ctof_p, float dt_ctof_pip, float momentum);
void Fill_deltat_after_mmsq_cut_Xpip(int pid, int charge, float dt_p, float dt_pip, float momentum);
void Fill_deltat_ctof_after_mmsq_cut_Xpip(int pid, int charge, float dt_ctof_p, float dt_ctof_pip, float momentum);

void Fill_theta_P_inv_mass(float W, float q2, float inv_mass, float theta, float wt);
void Fill_theta_pim_inv_mass(float W, float q2, float inv_mass, float theta, float wt);
void Fill_theta_pip_inv_mass(float W, float q2, float inv_mass, float theta, float wt);

void Fill_theta_P_lab_inv_mass(float W, float q2, float inv_mass, float theta, float wt);
void Fill_theta_pim_lab_inv_mass(float W, float q2, float inv_mass, float theta, float wt);
void Fill_theta_pip_lab_inv_mass(float W, float q2, float inv_mass, float theta, float wt);

void Fill_theta_P_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt);
void Fill_theta_pim_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt);
void Fill_theta_pip_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt);

void Fill_theta_P_lab_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt);
void Fill_theta_pim_lab_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt);
void Fill_theta_pip_lab_inv_mass_thrown(float W, float q2, float inv_mass, float theta, float wt);

void Fill_lu_dist(float li);
void Fill_lv_dist(float li);
void Fill_lw_dist(float li);

void Fill_vertex_vz(float vz);
void Fill_pid_size_fn(int pid_size, int pid);
void Fill_mom_ele(float mom, float mom_z);
void Fill_mom_p(float mom, float mom_z);
void Fill_mom_pip(float mom, float mom_z);
void Fill_mom_pim(float mom, float mom_z);
void Fill_mom_diff_e_pim(float mom_diff, float mom_diff_x, float mom_diff_y, float mom_diff_z);
void Fill_mom_diff_p_pip(float mom_diff, float mom_diff_x, float mom_diff_y, float mom_diff_z);
void Write_deltat();

// EC Sampling Fraction
void Fill_EC_sampling_fraction(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_1(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_2(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_3(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_4(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_5(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_6(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_7(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_8(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_9(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_10(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_11(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_12(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_13(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_14(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_15(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_16(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_17(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_18(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_19(float momentum, float sf, int sec_number, float weight);
void Fill_1d_sampling_fraction_20(float momentum, float sf, int sec_number, float weight);
void Fill_PCAL_VS_ECAL(float pcal, float ecal, int sec_number, float weight);
void Fill_hist_PCAL_FID_CUT(float x_PCAL, float y_PCAL);
void Fill_hist_PCAL_without_FID_CUT(float x_PCAL, float y_PCAL);
void Fill_hist_DC_FID_CUT(float R1X, float R1Y, float R2X, float R2Y, float R3X, float R3Y);
void Fill_hist_DC_without_FID_CUT(float R1X, float R1Y, float R2X, float R2Y, float R3X, float R3Y);
void Write_EC();
};

#endif
