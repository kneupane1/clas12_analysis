/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include "TLorentzVector.h"
#include "constants.hpp"
#include "physics.hpp"

class Reaction {
private:
TLorentzVector *_beam;
TLorentzVector *_elec;
TLorentzVector *_target;
TLorentzVector *_prot;
TLorentzVector *_pip;
TLorentzVector *_pim;

TLorentzVector *_q_cm;
TLorentzVector *_p_mu_prime_cm;
TLorentzVector *_pip_mu_prime_cm;
TLorentzVector *_pim_mu_prime_cm;

bool _hasE;
bool _hasP;
bool _hasPip;
bool _hasPim;

float _MM;
float _MM2;
float _MM_wop;
float _MM2_wop;

float _W;
float _Q2;

float _W_ep;
float _W_2pi;
float _W_delta_pp;
float _W_delta_zero;
float _W_rho;

float _W_singlepip;
float _Q2_2pi;
float _beta;
float _gamma;

float _theta_gamma = std::nan("-99");
float _phi_gamma = std::nan("-99");
float _theta_prot = std::nan("-99");
float _phi_prot = std::nan("-99");
float _theta_pip = std::nan("-99");
float _phi_pip = std::nan("-99");
float _theta_pim = std::nan("-99");
float _phi_pim = std::nan("-99");
float _alpha_ppip_pipim = std::nan("-99");
float _alpha_pippim_pipf = std::nan("-99");
float _alpha_ppim_pipip = std::nan("-99");

public:
Reaction();
Reaction(TLorentzVector *beam);
~Reaction();

void SetElec(float px, float py, float pz, float mass);
void SetProton(float px, float py, float pz, float mass);
void SetPip(float px, float py, float pz, float mass, int pid);
void SetPim(float px, float py, float pz, float mass);
TLorentzVector e_mu_prime();   // maile thapeko
TLorentzVector p_mu_prime();
TLorentzVector pip_mu_prime();
TLorentzVector pim_mu_prime();
//  TLorentzVector kp_mu_prime();
// TLorentzVector km_mu_prime();
// TLorentzVector q_cm(); // maile thapeko
float q_3_();
TLorentzVector p_mu_prime_cm();
TLorentzVector pip_mu_prime_cm();
TLorentzVector pim_mu_prime_cm();
double theta_();

//  void boost_fn(/*TLorentzVector four_vect, TLorentzVector e_mu,
// TLorentzVector e_mu_prime);
void CalcMissMass();
void CalcMissMass_wop();

void AlphaCalc();
float MM();
float MM2();
float MM_wop();
float MM2_wop();
float W();
float Q2();

float alpha_ppip_pipim();

float beta();
float gamma_();

float W_ep();
float W_2pi();
float W_delta_pp();
float W_delta_zero();
float W_rho();

float W_singlepip();
float Q2_2pi();
bool elecWopEvent();
bool twoPionWopEvent();
bool WopPimEvent();
bool WopPipEvent();
bool elecProtEvent();
bool twoPionEvent();
bool ProtonPimEvent();
bool ProtonPipEvent();
};

#endif