/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "physics.hpp"
#include "reaction.hpp"

Reaction::Reaction() {
  _beam = new TLorentzVector();
  _target = new TLorentzVector(0.0, 0.0, 0.0, MASS_P);
  _elec = new TLorentzVector();
  _prot = new TLorentzVector();
  _pip = new TLorentzVector();
  _pim = new TLorentzVector();

  _hasE = false;
  _hasP = false;
  _hasPip = false;
  _hasPim = false;

  _MM = std::nan("-99");
  _MM2 = std::nan("-99");
  _MM_wop = std::nan("-99");
  _MM2_wop = std::nan("-99");

  _W = std::nan("-99");
  _Q2 = std::nan("-99");
}
Reaction::Reaction(TLorentzVector *beam) {
  _beam = beam;
  _target = new TLorentzVector(0.0, 0.0, 0.0, MASS_P);
  _elec = new TLorentzVector();
  _prot = new TLorentzVector();
  _pip = new TLorentzVector();
  _pim = new TLorentzVector();

  _hasE = false;
  _hasP = false;
  _hasPip = false;
  _hasPim = false;

  _MM = std::nan("-99");
  _MM2 = std::nan("-99");
  _MM_wop = std::nan("-99");
  _MM2_wop = std::nan("-99");

  _W = std::nan("-99");
  _Q2 = std::nan("-99");
}
Reaction::~Reaction() {
  // delete _beam;
  // delete _elec;
  // delete _prot;
  // delete _pip;
  // delete _pim;
}

void Reaction::SetElec(float px, float py, float pz, float mass) {
  _hasE = true;
  _elec->SetXYZM(px, py, pz, mass);

  // Can calculate W and Q2 here
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
}

void Reaction::SetProton(float px, float py, float pz, float mass) {
  _hasP = true;
  _prot->SetXYZM(px, py, pz, mass);
}
void Reaction::SetPip(float px, float py, float pz, float mass) {
  _hasPip = true;
  _pip->SetXYZM(px, py, pz, mass);
}
void Reaction::SetPim(float px, float py, float pz, float mass) {
  _hasPim = true;
  _pim->SetXYZM(px, py, pz, mass);
}
// TLorentzVector Reaction::p_mu_prime() { return *_prot; }

void Reaction::CalcMissMass() {
  TLorentzVector mm;
  if (elecProtEvent()) {
    mm = (*_beam - *_elec);
    mm += *_target;
    mm -= *_prot;
    _MM = mm.M();
    _MM2 = mm.M2();
  } else if (twoPionEvent()) {
    mm = (*_beam - *_elec);
    mm += *_target;
    mm -= *_prot;
    mm -= *_pip;
    mm -= *_pim;
    _MM = mm.M();
    _MM2 = mm.M2();
  } else if (ProtonPimEvent()) {
    mm = (*_beam - *_elec);
    mm += *_target;
    mm -= *_prot;
    mm -= *_pim;
    _MM = mm.M();
    _MM2 = mm.M2();
  } else if (ProtonPipEvent()) {
    mm = (*_beam - *_elec);
    mm += *_target;
    mm -= *_prot;
    mm -= *_pip;
    _MM = mm.M();
    _MM2 = mm.M2();
  }
}
void Reaction::CalcMissMass_wop() {
  TLorentzVector mm_1;
  if (elecWopEvent()) {
    mm_1 = (*_beam - *_elec);
    mm_1 += *_target;
    // mm_1 -= *_prot;
    _MM_wop = mm_1.M();
    _MM2_wop = mm_1.M2();
  } else if (twoPionWopEvent()) {
    mm_1 = (*_beam - *_elec);
    mm_1 += *_target;
    // mm_1 -= *_prot;
    mm_1 -= *_pip;
    mm_1 -= *_pim;
    _MM_wop = mm_1.M();
    _MM2_wop = mm_1.M2();
  } else if (WopPimEvent()) {
    mm_1 = (*_beam - *_elec);
    mm_1 += *_target;
    // mm_1 -= *_prot;
    mm_1 -= *_pim;
    _MM_wop = mm_1.M();
    _MM2_wop = mm_1.M2();
  } else if (WopPipEvent()) {
    mm_1 = (*_beam - *_elec);
    mm_1 += *_target;
    // mm_1 -= *_prot;
    mm_1 -= *_pip;
    _MM_wop = mm_1.M();
    _MM2_wop = mm_1.M2();
  }
}
TLorentzVector Reaction::e_mu_prime() { return *_elec; } // maile thapeko
TLorentzVector Reaction::p_mu_prime() { return *_prot; }
TLorentzVector Reaction::pip_mu_prime() { return *_pip; }
TLorentzVector Reaction::pim_mu_prime() { return *_pim; }

float Reaction::MM() { return _MM; }
float Reaction::MM2() { return _MM2; }
float Reaction::MM_wop() { return _MM_wop; }
float Reaction::MM2_wop() { return _MM2_wop; }

float Reaction::W() { return _W; }
float Reaction::Q2() { return _Q2; }

bool Reaction::elecProtEvent() {
  return (_hasE && _hasP && !_hasPip && !_hasPim);
}
bool Reaction::twoPionEvent() { return (_hasE && _hasP && _hasPip && _hasPim); }
bool Reaction::ProtonPimEvent() {
  return (_hasE && _hasP && _hasPim && !_hasPip);
}
bool Reaction::ProtonPipEvent() {
  return (_hasE && _hasP && _hasPip && !_hasPim);
}
bool Reaction::elecWopEvent() {
  return (_hasE /*&& _hasP*/ && !_hasPip && !_hasPim);
}
bool Reaction::twoPionWopEvent() {
  return (_hasE /*&& _hasP*/ && _hasPip && _hasPim);
}
bool Reaction::WopPimEvent() {
  return (_hasE /*&& _hasP*/ && _hasPim && !_hasPip);
}
bool Reaction::WopPipEvent() {
  return (_hasE /*&& _hasP*/ && _hasPip && !_hasPim);
}
