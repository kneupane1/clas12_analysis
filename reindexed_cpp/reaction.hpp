

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
TLorentzVector *_neutron;
TLorentzVector *_other;
TLorentzVector *_elec_thrown;
TLorentzVector *_prot_thrown;
TLorentzVector *_pip_thrown;
TLorentzVector *_pim_thrown;

TLorentzVector *_q_cm;
TLorentzVector *_p_mu_prime_cm;
TLorentzVector *_pip_mu_prime_cm;
TLorentzVector *_pim_mu_prime_cm;

TLorentzVector *_q_cm_thrown;
TLorentzVector *_p_mu_prime_cm_thrown;
TLorentzVector *_pip_mu_prime_cm_thrown;
TLorentzVector *_pim_mu_prime_cm_thrown;

bool _hasE = false;
bool _hasP = false;
bool _hasPip = false;
bool _hasPim = false;
bool _hasOther = false;
bool _hasNeutron = false;

short _numProt = 0;
short _numPip = 0;
short _numPim = 0;
short _numPos = 0;
short _numNeg = 0;
short _numNeutral = 0;
short _numOther = 0;

float _MM;
float _MM2;
float _MM_wop;
float _MM2_wop;

float _W;
float _Q2;

float _W_thrown;
float _Q2_thrown;

float _W_ep;
float _W_2pi;
float _W_delta_pp;
float _W_delta_zero;
float _W_rho;

float _W_2pi_thrown;
float _W_delta_pp_thrown;
float _W_delta_zero_thrown;
float _W_rho_thrown;

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

float _theta_gamma_thrown = std::nan("-99");
float _phi_gamma_thrown = std::nan("-99");
float _theta_prot_thrown = std::nan("-99");
float _phi_prot_thrown = std::nan("-99");
float _theta_pip_thrown = std::nan("-99");
float _phi_pip_thrown = std::nan("-99");
float _theta_pim_thrown = std::nan("-99");
float _phi_pim_thrown = std::nan("-99");
float _alpha_ppip_pipim_thrown = std::nan("-99");
float _alpha_pippim_pipf_thrown = std::nan("-99");
float _alpha_ppim_pipip_thrown = std::nan("-99");
public:
Reaction();
Reaction(TLorentzVector *beam);
~Reaction();

void SetElec(float px, float py, float pz, float mass);
void SetProton(float px, float py, float pz, float mass);
void SetPip(float px, float py, float pz, float mass);
void SetPim(float px, float py, float pz, float mass);
void SetOther(float px, float py, float pz, float mass, int pid);

void SetElec_thrown(float px, float py, float pz, float mass);
void SetProt_thrown(float px, float py, float pz, float mass);
void SetPip_thrown(float px, float py, float pz, float mass);
void SetPim_thrown(float px, float py, float pz, float mass);
TLorentzVector e_mu_prime();    // maile thapeko
TLorentzVector p_mu_prime();
TLorentzVector pip_mu_prime();
TLorentzVector pim_mu_prime();

TLorentzVector e_mu_prime_thrown();    // maile thapeko
TLorentzVector p_mu_prime_thrown();
TLorentzVector pip_mu_prime_thrown();
TLorentzVector pim_mu_prime_thrown();

//  TLorentzVector kp_mu_prime();
// TLorentzVector km_mu_prime();
// TLorentzVector q_cm(); // maile thapeko
float q_3_();
TLorentzVector p_mu_prime_cm();
TLorentzVector pip_mu_prime_cm();
TLorentzVector pim_mu_prime_cm();
//float theta_();
TLorentzVector p_mu_prime_cm_thrown();
TLorentzVector pip_mu_prime_cm_thrown();
TLorentzVector pim_mu_prime_cm_thrown();
//  void boost_fn(/*TLorentzVector four_vect, TLorentzVector e_mu,
// TLorentzVector e_mu_prime);
void CalcMissMass();
void CalcMissMass_wop();
void thrownCalc();

void AlphaCalc();
void AlphaCalc_thrown();
float MM();
float MM2();
float MM_wop();
float MM2_wop();
float W();
float Q2();

float W_thrown();
float Q2_thrown();

float alpha_ppip_pipim();
float alpha_pippim_pipf();
float alpha_ppim_pipip();

float alpha_ppip_pipim_thrown();
float alpha_pippim_pipf_thrown();
float alpha_ppim_pipip_thrown();

float beta();
float gamma_();

float W_ep();
float W_2pi();
float W_delta_pp();
float W_delta_zero();
float W_rho();

float W_2pi_thrown();
float W_delta_pp_thrown();
float W_delta_zero_thrown();
float W_rho_thrown();

float W_singlepip();
float Q2_2pi();
// bool elecWopEvent();
// bool twoPionWopEvent();
// bool WopPimEvent();
// bool WopPipEvent();
// bool elecProtEvent();
// bool twoPionEvent();
// bool ProtonPimEvent();
// bool ProtonPipEvent();
bool elecProtEvent() {
        return ((_numProt == 1) && (_hasE && _hasP && !_hasPip && !_hasPim));
}
bool twoPionEvent() {
        return ((/*_numProt == 1 &&*/ _numPip == 1 && _numPim == 1) && (_hasE && _hasP && _hasPip && _hasPim  && !_hasOther));
}
bool ProtonPimEvent() {
        return ((_numProt == 1 && _numPim == 1) && (_hasE && _hasP && _hasPim && !_hasPip));
}
bool ProtonPipEvent() {
        return ((_numProt == 1 && _numPip == 1) && (_hasE && _hasP && _hasPip && !_hasPim));
}
bool elecWopEvent() {
        return (_hasE /*&& _hasP*/ && !_hasPip && !_hasPim);
}
bool twoPionWopEvent() {
        return ((_numProt == 1 && _numPim == 1) && (_hasE && _hasPip && _hasPim));
}
bool WopPimEvent() {
        return ((_numPim == 1) && (_hasE /*&& _hasP*/ && _hasPim && !_hasPip));
}
bool WopPipEvent() {
        return ((_numPip == 1) && (_hasE && !_hasP && _hasPip && !_hasPim));
}

// extra
bool elecPimEvent() {
        return ((_numPim == 1) && (_hasE && !_hasP && !_hasPip && _hasPim));
}

};

#endif
