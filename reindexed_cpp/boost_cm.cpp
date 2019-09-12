#include "boost.hpp"
#include "physics.hpp"
using namespace std;
Boost_cm::Boost_cm() {
  _beam = new TLorentzVector();
  _target = new TLorentzVector(0.0, 0.0, 0.0, MASS_P);
  _elec = new TLorentzVector();
  _prot = new TLorentzVector();
  _pip = new TLorentzVector();
  _pim = new TLorentzVector();

  _q_cm = new TLorentzVector();
  _p_mu_prime_cm = new TLorentzVector();
  _pip_mu_prime_cm = new TLorentzVector();
  _pim_mu_prime_cm = new TLorentzVector();

  _hasE = false;
  _hasP = false;
  _hasPip = false;
  _hasPim = false;
}
Boost_cm::Boost_cm(TLorentzVector *beam) {
  _beam = beam;
  _target = new TLorentzVector(0.0, 0.0, 0.0, MASS_P);
  _elec = new TLorentzVector();
  _prot = new TLorentzVector();
  _pip = new TLorentzVector();
  _pim = new TLorentzVector();

  _q_cm = new TLorentzVector();
  _p_mu_prime_cm = new TLorentzVector();
  _pip_mu_prime_cm = new TLorentzVector();
  _pim_mu_prime_cm = new TLorentzVector();

  _hasE = false;
  _hasP = false;
  _hasPip = false;
  _hasPim = false;
}
Boost_cm::~Boost_cm() {
  // delete _beam;
  delete _elec;
  delete _prot;
  delete _pip;
  delete _pim;
  delete _target;

  delete _q_cm;
  delete _p_mu_prime_cm;
  delete _pip_mu_prime_cm;
  delete _pim_mu_prime_cm;
}
bool Reaction::twoPionEvent() { return (_hasE && _hasP && _hasPip && _hasPim); }

float Reaction::q_cm() {
  if (twoPionEvent())
    return (physics::boost_((*_beam - *_elec), *_beam, *_elec));
}
double Reaction::theta_() {
  return physics::theta_fn(*_pip_mu_prime_cm, *_beam, *_elec);
}
TLorentzVector Reaction::p_mu_prime_cm() {
  if (twoPionEvent())
    return (physics::boost_(*_prot, *_beam, *_elec)); // /*_p_mu_prime_cm; /*/
}
TLorentzVector Reaction::pip_mu_prime_cm() {
  if (twoPionEvent())
    return *_pip_mu_prime_cm; // (physics::boost_(*_pip, *_beam, *_elec));
}
TLorentzVector Reaction::pim_mu_prime_cm() {
  if (twoPionEvent())
    return *_pim_mu_prime_cm; //(physics::boost_(*_pim, *_beam, *_elec));
}
void PhiAng(float phi) {
  if (phi > 0)
    return phi = phi;
  else if (phi < 0)
    return phi = phi + 360;
}
if (pim_mu_prime_cm().Phi() > 0)
  phi_PIm_cm = (180. / PI) * pim_mu_prime_cm().Phi();
if (pim_mu_prime_cm().Phi() < 0)
  phi_PIm_cm = (180. / PI) * (pim_mu_prime_cm().Phi() + 2 * PI);

if (pip_mu_prime_cm().Phi() > 0)
  phi_PIp_cm = (180. / PI) * pip_mu_prime_cm().Phi();
if (pip_mu_prime_cm().Phi() < 0)
  phi_PIp_cm = (180. / PI) * (pip_mu_prime_cm().Phi() + 2 * PI);

if (p_mu_prime_cm().Phi() > 0)
  phi_P_cm = (180. / PI) * p_mu_prime_cm().Phi();
if (p_mu_prime_cm().Phi() < 0)
  phi_P_cm = (180. / PI) * (p_mu_prime_cm().Phi() + 2 * PI);
}
/// 1
void AlphaCalc() {
  a_gamma =
      sqrt(1. / (1 - pow((pim_mu_prime_cm().Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(pim_mu_prime_cm().Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * pim_mu_prime_cm().Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((pim_mu_prime_cm().Vect().Unit() *
                               pip_mu_prime_cm().Vect().Unit()),
                              2)));
  b_beta =
      -(pim_mu_prime_cm().Vect().Unit() * pip_mu_prime_cm().Vect().Unit()) *
      a_beta;
  Vect3_beta = a_beta * pip_mu_prime_cm().Vect().Unit() +
               b_beta * pim_mu_prime_cm().Vect().Unit();

  alpha_PPIp_piPIm = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
  if (Vect3_gamma.Cross(Vect3_beta) * pim_mu_prime_cm().Vect() < 0)
    alpha_PPIp_piPIm = 360. - alpha_PPIp_piPIm;
}
