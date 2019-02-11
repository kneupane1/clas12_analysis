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
  float _W_singlepip;
  float _Q2_2pi;

public:
  Reaction();
  Reaction(TLorentzVector *beam);
  ~Reaction();

  void SetElec(float px, float py, float pz, float mass);
  void SetProton(float px, float py, float pz, float mass);
  void SetPip(float px, float py, float pz, float mass);
  void SetPim(float px, float py, float pz, float mass);
  TLorentzVector e_mu_prime(); // maile thapeko
  TLorentzVector p_mu_prime();
  TLorentzVector pip_mu_prime();
  TLorentzVector pim_mu_prime();
  //  TLorentzVector kp_mu_prime();
  // TLorentzVector km_mu_prime();

  void CalcMissMass();
  void CalcMissMass_wop();

  float MM();
  float MM2();
  float MM_wop();
  float MM2_wop();
  float W();
  float Q2();
  float W_ep();
  float W_2pi();
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
