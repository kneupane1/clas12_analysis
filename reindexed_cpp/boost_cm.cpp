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

  _theta_gamma = std::nan("-99");
  _phi_gamma = std::nan("-99");
  _theta_prot = std::nan("-99");
  _phi_prot = std::nan("-99");
  _theta_pip = std::nan("-99");
  _phi_pip = std::nan("-99");
  _theta_pim = std::nan("-99");
  _phi_pim = std::nan("-99");
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

  _theta_gamma = std::nan("-99");
  _phi_gamma = std::nan("-99");
  _theta_prot = std::nan("-99");
  _phi_prot = std::nan("-99");
  _theta_pip = std::nan("-99");
  _phi_pip = std::nan("-99");
  _theta_pim = std::nan("-99");
  _phi_pim = std::nan("-99");

  _alpha_pip_pim = std::nan("-99");
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
bool Boost_cm::twoPionEvent() { return (_hasE && _hasP && _hasPip && _hasPim); }

float Boost_cm::AlphaCalc() {
  TLorentzVector q_cm_;
  double theta_gamma;
  double phi_gamma;
  double theta_prot;
  double phi_prot;
  double theta_pip;
  double phi_pip;
  double theta_pim;
  double phi_pim;

  //  Float_t m_proton, m_pip, beta;
  Float_t a_gamma, b_gamma, a_beta, b_beta;
  TVector3 Vect3_gamma, Vect3_beta, V3_anti_z(0, 0, -1);
  float alpha_PPIp_piPIm;
  if (twoPionEvent()) {
    *_q_cm = (physics::boost_((*_beam - *_elec), *_beam, *_elec));
    *_p_mu_prime_cm = physics::boost_(*_prot, *_beam, *_elec);
    *_pip_mu_prime_cm = physics::boost_(*_pip, *_beam, *_elec);
    *_pim_mu_prime_cm = physics::boost_(*_pim, *_beam, *_elec);
    theta_gamma = *_q_cm.Theta() * 180 / PI;
    theta_prot = *_p_mu_prime_cm.Theta() * 180 / PI;
    theta_pip = *_pip_mu_prime_cm_.Theta() * 180 / PI;
    theta_pim = *_pim_mu_prime_cm_.Theta() * 180 / PI;

    if (*_q_cm.Phi() > 0)
      phi_gamma = *_q_cm_.Phi() * 180 / PI;
    else if (*_q_cm.Phi() < 0)
      phi_gamma = (*_q_cm.Phi() + 2 * Pi) * 180 / PI;

    if (*_p_mu_prime_cm.Phi() > 0)
      phi_gamma = *_p_mu_prime_cm_.Phi() * 180 / PI;
    else if (*_p_mu_prime_cm.Phi() < 0)
      phi_gamma = (*_p_mu_prime_cm.Phi() + 2 * Pi) * 180 / PI;
    if (*_pip_mu_prime_cm.Phi() > 0)
      phi_gamma = *_pip_mu_prime_cm_.Phi() * 180 / PI;
    else if (*_pip_mu_prime_cm.Phi() < 0)
      phi_gamma = (*_pip_mu_prime_cm.Phi() + 2 * Pi) * 180 / PI;
    if (*_pim_mu_prime_cm.Phi() > 0)
      phi_gamma = *_pim_mu_prime_cm_.Phi() * 180 / PI;
    else if (*_pim_mu_prime_cm.Phi() < 0)
      phi_gamma = (*_pim_mu_prime_cm.Phi() + 2 * Pi) * 180 / PI;

    a_gamma =
        sqrt(1. / (1 - pow((*_pim_mu_prime_cm.Vect().Unit() * V3_anti_z), 2)));
    b_gamma = -(*_pim_mu_prime_cm.Vect().Unit() * V3_anti_z) * a_gamma;
    Vect3_gamma =
        a_gamma * V3_anti_z + b_gamma * *_pim_mu_prime_cm.Vect().Unit();

    a_beta = sqrt(1. / (1 - pow((*_pim_mu_prime_cm.Vect().Unit() *
                                 *_pip_mu_prime_cm.Vect().Unit()),
                                2)));
    b_beta =
        -(*_pim_mu_prime_cm.Vect().Unit() * pip_mu_prime_cm().Vect().Unit()) *
        a_beta;
    Vect3_beta = a_beta * *_pip_mu_prime_cm.Vect().Unit() +
                 b_beta * *_pim_mu_prime_cm.Vect().Unit();

    alpha_PPIp_piPIm = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
    if (Vect3_gamma.Cross(Vect3_beta) * *_pim_mu_prime_cm.Vect() < 0)
      alpha_PPIp_piPIm = 360. - alpha_PPIp_piPIm;

    _theta_gamma = theta_gamma;
    _phi_gamma = phi_gamma;
    _theta_prot = theta_prot;
    _phi_prot = phi_prot;
    _theta_pip = theta_pip;
    _phi_pip = phi_pip;
    _theta_pim = theta_pim;
    _phi_pim = phi_pim;

    return alpha_PPIp_piPIm;
  }
}
// double Reaction::theta_() {
//   return physics::theta_fn(*_pip_mu_prime_cm, *_beam, *_elec);
// }
// TLorentzVector Reaction::p_mu_prime_cm() {
//   if (twoPionEvent())
//     return (physics::boost_(*_prot, *_beam, *_elec)); // /*_p_mu_prime_cm;
//     /*/
// }
// TLorentzVector Reaction::pip_mu_prime_cm() {
//   if (twoPionEvent())
//     return *_pip_mu_prime_cm; // (physics::boost_(*_pip, *_beam, *_elec));
// }
// TLorentzVector Reaction::pim_mu_prime_cm() {
//   if (twoPionEvent())
//     return *_pim_mu_prime_cm; //(physics::boost_(*_pim, *_beam, *_elec));
// }
// void PhiAng(float phi) {
//   if (phi > 0)
//     return phi = phi;
//   else if (phi < 0)
//     return phi = phi + 360;
// }
//
// if (pim_mu_prime_cm().Phi() > 0)
//   phi_PIm_cm = (180. / PI) * pim_mu_prime_cm().Phi();
// if (pim_mu_prime_cm().Phi() < 0)
//   phi_PIm_cm = (180. / PI) * (pim_mu_prime_cm().Phi() + 2 * PI);
//
// if (pip_mu_prime_cm().Phi() > 0)
//   phi_PIp_cm = (180. / PI) * pip_mu_prime_cm().Phi();
// if (pip_mu_prime_cm().Phi() < 0)
//   phi_PIp_cm = (180. / PI) * (pip_mu_prime_cm().Phi() + 2 * PI);
//
// if (p_mu_prime_cm().Phi() > 0)
//   phi_P_cm = (180. / PI) * p_mu_prime_cm().Phi();
// if (p_mu_prime_cm().Phi() < 0)
//   phi_P_cm = (180. / PI) * (p_mu_prime_cm().Phi() + 2 * PI);
// }
// /// 1
// void AlphaCalc() {
//   a_gamma =
//       sqrt(1. / (1 - pow((pim_mu_prime_cm().Vect().Unit() * V3_anti_z), 2)));
//   b_gamma = -(pim_mu_prime_cm().Vect().Unit() * V3_anti_z) * a_gamma;
//   Vect3_gamma = a_gamma * V3_anti_z + b_gamma *
//   pim_mu_prime_cm().Vect().Unit();
//
//   a_beta = sqrt(1. / (1 - pow((pim_mu_prime_cm().Vect().Unit() *
//                                pip_mu_prime_cm().Vect().Unit()),
//                               2)));
//   b_beta =
//       -(pim_mu_prime_cm().Vect().Unit() * pip_mu_prime_cm().Vect().Unit()) *
//       a_beta;
//   Vect3_beta = a_beta * pip_mu_prime_cm().Vect().Unit() +
//                b_beta * pim_mu_prime_cm().Vect().Unit();
//
//   alpha_PPIp_piPIm = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
//   if (Vect3_gamma.Cross(Vect3_beta) * pim_mu_prime_cm().Vect() < 0)
//     alpha_PPIp_piPIm = 360. - alpha_PPIp_piPIm;
// }
