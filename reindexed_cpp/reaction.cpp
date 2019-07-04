
#include "physics.hpp"
#include "reaction.hpp"
using namespace std;
Reaction::Reaction() {
  _beam = new TLorentzVector();
  _target = new TLorentzVector(0.0, 0.0, 0.0, MASS_P);
  _elec = new TLorentzVector();
  _prot = new TLorentzVector();
  _pip = new TLorentzVector();
  _pim = new TLorentzVector();
  _neutron = new TLorentzVector();
  _other = new TLorentzVector();

  _q_cm = new TLorentzVector();
  _p_mu_prime_cm = new TLorentzVector();
  _pip_mu_prime_cm = new TLorentzVector();
  _pim_mu_prime_cm = new TLorentzVector();

  _MM = std::nan("-99");
  _MM2 = std::nan("-99");
  _MM_wop = std::nan("-99");
  _MM2_wop = std::nan("-99");

  _W = std::nan("-99");
  _Q2 = std::nan("-99");
  _W_ep = std::nan("-99");
  _W_2pi = std::nan("-99");
  _W_delta_pp = std::nan("-99");
  _W_delta_zero = std::nan("-99");
  _W_rho = std::nan("-99");

  _W_singlepip = std::nan("-99");

  _Q2_2pi = std::nan("-99");

  _beta = std::nan("-99");
  _gamma = std::nan("-99");

  _theta_gamma = std::nan("-99");
  _phi_gamma = std::nan("-99");
  _theta_prot = std::nan("-99");
  _phi_prot = std::nan("-99");
  _theta_pip = std::nan("-99");
  _phi_pip = std::nan("-99");
  _theta_pim = std::nan("-99");
  _phi_pim = std::nan("-99");
  _alpha_ppip_pipim = std::nan("-99");
  _alpha_pippim_pipf = std::nan("-99");
  _alpha_ppim_pipip = std::nan("-99");
}
Reaction::Reaction(TLorentzVector *beam) {
  _beam = beam;
  _target = new TLorentzVector(0.0, 0.0, 0.0, MASS_P);
  _elec = new TLorentzVector();
  _prot = new TLorentzVector();
  _pip = new TLorentzVector();
  _pim = new TLorentzVector();
  _neutron = new TLorentzVector();
  _other = new TLorentzVector();

  _q_cm = new TLorentzVector();
  _p_mu_prime_cm = new TLorentzVector();
  _pip_mu_prime_cm = new TLorentzVector();
  _pim_mu_prime_cm = new TLorentzVector();

  // _hasE = false;
  // _hasP = false;
  // _hasPip = false;
  // _hasPim = false;
  // _hasNeutron = false;
  // _hasOther = false;

  _MM = std::nan("-99");
  _MM2 = std::nan("-99");
  _MM_wop = std::nan("-99");
  _MM2_wop = std::nan("-99");

  _W = std::nan("-99");
  _Q2 = std::nan("-99");

  _W_ep = std::nan("-99");
  _W_2pi = std::nan("-99");
  _W_delta_pp = std::nan("-99");
  _W_delta_zero = std::nan("-99");
  _W_rho = std::nan("-99");

  _W_singlepip = std::nan("-99");
  _Q2_2pi = std::nan("-99");

  _beta = std::nan("-99");
  _gamma = std::nan("-99");

  _theta_gamma = std::nan("-99");
  _phi_gamma = std::nan("-99");
  _theta_prot = std::nan("-99");
  _phi_prot = std::nan("-99");
  _theta_pip = std::nan("-99");
  _phi_pip = std::nan("-99");
  _theta_pim = std::nan("-99");
  _phi_pim = std::nan("-99");
  _alpha_ppip_pipim = std::nan("-99");
  _alpha_pippim_pipf = std::nan("-99");
  _alpha_ppim_pipip = std::nan("-99");
}
Reaction::~Reaction() {
  // delete _beam;
  delete _elec;
  delete _prot;
  delete _pip;
  delete _pim;
  delete _target;
  delete _neutron;
  delete _other;

  delete _q_cm;
  delete _p_mu_prime_cm;
  delete _pip_mu_prime_cm;
  delete _pim_mu_prime_cm;
}

void Reaction::SetElec(float px, float py, float pz, float mass) {
  _hasE = true;
  _elec->SetXYZM(px, py, pz, mass);

  _q_cm->SetXYZM(0.0, 0.0, 0.0, 0.0);
  _p_mu_prime_cm->SetXYZM(0.0, 0.0, 0.0, 0.0);
  _pip_mu_prime_cm->SetXYZM(0.0, 0.0, 0.0, 0.0);
  _pim_mu_prime_cm->SetXYZM(0.0, 0.0, 0.0, 0.0);

  // Can calculate W and Q2 here
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
  _beta = physics::beta_boost(*_beam, *_elec);
  _gamma = physics::gamma(*_beam, *_elec);
}

void Reaction::SetProton(float px, float py, float pz, float mass) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _prot->SetXYZM(px, py, pz, mass);
}
void Reaction::SetPip(float px, float py, float pz, float mass) {
  {
    _numPip++;
    _numPos++;
    _hasPip = true;
    _pip->SetXYZM(px, py, pz, mass);
  }
}
void Reaction::SetPim(float px, float py, float pz, float mass) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim->SetXYZM(px, py, pz, mass);
}
void Reaction::SetOther(float px, float py, float pz, float mass, int pid) {
  if (pid == NEUTRON) {
    _hasNeutron = true;
    _numNeutral++;
    _neutron->SetXYZM(px, py, pz, mass);
  } else {
    _hasOther = true;
    _numOther++;
    _other->SetXYZM(px, py, pz, mass);
  }
}
void Reaction::CalcMissMass() {
  TLorentzVector mm;
  TLorentzVector m_dpp;
  TLorentzVector m_d0;
  TLorentzVector m_rho;
  TLorentzVector w_2pi;

  mm = (*_beam - *_elec + *_target);
  if (elecProtEvent()) {
    mm -= *_prot;
    _MM = mm.M();
    _MM2 = mm.M2();
    _W_ep = physics::W_calc(*_beam, *_elec);

  } else if (twoPionEvent()) {
    mm -= *_prot;
    mm -= *_pip;
    mm -= *_pim;
    _MM = mm.M();
    _MM2 = mm.M2();
    //  if (twoPionEvent()) {
    w_2pi = (*_prot);
    w_2pi += *_pip;
    w_2pi += *_pim;
    _W_2pi = w_2pi.M();  // invariant mass pf P' pip pim

    //  _W_2pi = physics::W_calc(*_beam, *_elec);
    _Q2_2pi = physics::Q2_calc(*_beam, *_elec);
    m_dpp = (*_prot);
    m_dpp += *_pip;
    _W_delta_pp = m_dpp.M();  // invariant mass of P' pip
    m_d0 = (*_prot);
    m_d0 += *_pim;
    _W_delta_zero = m_d0.M();  // invariant mass of P' pim
    m_rho = (*_pip);
    m_rho += *_pim;
    _W_rho = m_rho.M();  // invariant mass of pip pim
    *_pip_mu_prime_cm = physics::boost_(*_pip, *_beam, *_elec);
    // *_q_cm = physics::boost_((*_beam - *_elec), *_beam, *_elec);
    *_p_mu_prime_cm = physics::boost_(*_prot, *_beam, *_elec);
    *_pim_mu_prime_cm = physics::boost_(*_pim, *_beam, *_elec);
    //}
  } else if (ProtonPimEvent()) {
    mm -= *_prot;
    mm -= *_pim;
    _MM = mm.M();
    _MM2 = mm.M2();
  } else if (ProtonPipEvent()) {
    mm -= *_prot;
    mm -= *_pip;
    _MM = mm.M();
    _MM2 = mm.M2();
  }
}
void Reaction::CalcMissMass_wop() {
  TLorentzVector mm_1;
  mm_1 = (*_beam - *_elec + *_target);
  if (elecWopEvent()) {
    _MM_wop = mm_1.M();
    _MM2_wop = mm_1.M2();
  } else if (twoPionWopEvent()) {
    mm_1 -= *_pip;
    mm_1 -= *_pim;
    _MM_wop = mm_1.M();
    _MM2_wop = mm_1.M2();
  } else if (WopPimEvent()) {
    mm_1 -= *_pim;
    _MM_wop = mm_1.M();
    _MM2_wop = mm_1.M2();
  } else if (WopPipEvent()) {
    mm_1 -= *_pip;
    _MM_wop = mm_1.M();
    _MM2_wop = mm_1.M2();
    _W_singlepip = physics::W_calc(*_beam, *_elec);
  }
}
TLorentzVector Reaction::e_mu_prime() { return *_elec; }  // maile thapeko
TLorentzVector Reaction::p_mu_prime() { return *_prot; }
TLorentzVector Reaction::pip_mu_prime() { return *_pip; }
TLorentzVector Reaction::pim_mu_prime() { return *_pim; }

// TLorentzVector Reaction::q_cm() {
//  if (twoPionEvent())
float Reaction::q_3_() {
  return physics::q_3(*_beam, *_elec);  // *_q_cm;  (physics::boost_((*_beam -
                                        // *_elec), *_beam, *_elec));
}
// float Reaction::theta_() {
//         return physics::theta_fn(*_pip_mu_prime_cm, *_beam, *_elec);
// }

TLorentzVector Reaction::p_mu_prime_cm() {
  // if (twoPionEvent())
  return *_p_mu_prime_cm; /*
                             (physics::boost_(*_prot, *_beam, *_elec));*/
}
TLorentzVector Reaction::pip_mu_prime_cm() {
  // if (twoPionEvent())
  return *_pip_mu_prime_cm;  // (physics::boost_(*_pip, *_beam, *_elec));
}
TLorentzVector Reaction::pim_mu_prime_cm() {
  // if (twoPionEvent())
  return *_pim_mu_prime_cm;  //(physics::boost_(*_pim, *_beam, *_elec));
}
// void Reaction::AlphaCalc() {
//         TLorentzVector _q_cm_;
//         TLorentzVector _p_mu_prime_cm_;
//         TLorentzVector _pip_mu_prime_cm_;
//         TLorentzVector _pim_mu_prime_cm_;
//
//         float theta_gamma;
//         float phi_gamma;
//         float theta_prot;
//         float phi_prot;
//         float theta_pip;
//         float phi_pip;
//         float theta_pim;
//         float phi_pim;
//
//         //  Float_t m_proton, m_pip, beta;
//         Float_t a_gamma, b_gamma, a_beta, b_beta;
//         TVector3 Vect3_gamma, Vect3_beta, V3_anti_z(0, 0, -1);
//         float alpha_PPIp_piPIm;
//         float alpha_PIpPIm_pipf;
//         float alpha_PPIm_piPIp;
//
//         if (twoPionEvent()) {
//                 *_q_cm = (physics::boost_((*_beam - *_elec), *_beam, *_elec));
//                 // *_p_mu_prime_cm = physics::boost_(*_prot, *_beam, *_elec);
//                 // *_pip_mu_prime_cm = physics::boost_(*_pip, *_beam, *_elec);
//                 // *_pim_mu_prime_cm = physics::boost_(*_pim, *_beam, *_elec);
//                 _q_cm_ = (*_q_cm);
//                 _p_mu_prime_cm_ = (*_p_mu_prime_cm);
//                 _pip_mu_prime_cm_ = (*_pip_mu_prime_cm);
//                 _pim_mu_prime_cm_ = (*_pim_mu_prime_cm);
//
//                 theta_gamma = _q_cm_.Theta() * 180 / PI;
//                 theta_prot = _p_mu_prime_cm_.Theta() * 180 / PI;
//                 theta_pip = _pip_mu_prime_cm_.Theta() * 180 / PI;
//                 theta_pim = _pim_mu_prime_cm_.Theta() * 180 / PI;
//                 if (_q_cm_.Phi() > 0)
//                         phi_gamma = _q_cm_.Phi() * 180 / PI;
//                 else if (_q_cm_.Phi() < 0)
//                         phi_gamma = (_q_cm_.Phi() + 2 * PI) * 180 / PI;
//
//                 if (_p_mu_prime_cm_.Phi() > 0)
//                         phi_gamma = _p_mu_prime_cm_.Phi() * 180 / PI;
//                 else if (_p_mu_prime_cm_.Phi() < 0)
//                         phi_gamma = (_p_mu_prime_cm_.Phi() + 2 * PI) * 180 / PI;
//                 if (_pip_mu_prime_cm_.Phi() > 0)
//                         phi_gamma = _pip_mu_prime_cm_.Phi() * 180 / PI;
//                 else if (_pip_mu_prime_cm_.Phi() < 0)
//                         phi_gamma = (_pip_mu_prime_cm_.Phi() + 2 * PI) * 180 / PI;
//                 if (_pim_mu_prime_cm_.Phi() > 0)
//                         phi_gamma = _pim_mu_prime_cm_.Phi() * 180 / PI;
//                 else if (_pim_mu_prime_cm_.Phi() < 0)
//                         phi_gamma = (_pim_mu_prime_cm_.Phi() + 2 * PI) * 180 / PI;
//                 // 1 this one is used for α[π−]
//                 a_gamma = sqrt(1. / (1 - pow((_pim_mu_prime_cm_.Vect().Unit() * V3_anti_z),
//                                              2))); // V3_anti_z(0,0,-1);
//                 b_gamma = -(_pim_mu_prime_cm_.Vect().Unit() * V3_anti_z) * a_gamma;
//                 Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _pim_mu_prime_cm_.Vect().Unit();
//
//                 a_beta = sqrt(1. / (1 - pow((_pim_mu_prime_cm_.Vect().Unit() * _pip_mu_prime_cm_.Vect().Unit()),
//                 2))); b_beta = -(_pim_mu_prime_cm_.Vect().Unit() * _pip_mu_prime_cm_.Vect().Unit()) * a_beta;
//                 Vect3_beta = a_beta * _pip_mu_prime_cm_.Vect().Unit() + b_beta * _pim_mu_prime_cm_.Vect().Unit();
//
//                 alpha_PPIp_piPIm = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
//                 if (Vect3_gamma.Cross(Vect3_beta) * _pim_mu_prime_cm_.Vect() < 0)
//                         alpha_PPIp_piPIm = 360. - alpha_PPIp_piPIm; //α[pπ+][p'π−]
//
//                 /// 2
//                 a_gamma = sqrt(1. / (1 - pow((_p_mu_prime_cm_.Vect().Unit() * V3_anti_z), 2)));
//                 b_gamma = -(_p_mu_prime_cm_.Vect().Unit() * V3_anti_z) * a_gamma;
//                 Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _p_mu_prime_cm_.Vect().Unit();
//
//                 a_beta = sqrt(1. / (1 - pow((_p_mu_prime_cm_.Vect().Unit() * _pip_mu_prime_cm_.Vect().Unit()), 2)));
//                 b_beta = -(_p_mu_prime_cm_.Vect().Unit() * _pip_mu_prime_cm_.Vect().Unit()) * a_beta;
//                 Vect3_beta = a_beta * _pip_mu_prime_cm_.Vect().Unit() + b_beta * _p_mu_prime_cm_.Vect().Unit();
//
//                 alpha_PIpPIm_pipf = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
//
//                 if (Vect3_gamma.Cross(Vect3_beta) * _p_mu_prime_cm_.Vect() < 0)
//                         alpha_PIpPIm_pipf = 360. - alpha_PIpPIm_pipf; //α[pp'][π+π−]
//
//                 /// 3
//                 a_gamma = sqrt(1. / (1 - pow((_pip_mu_prime_cm_.Vect().Unit() * V3_anti_z), 2)));
//                 b_gamma = -(_pip_mu_prime_cm_.Vect().Unit() * V3_anti_z) * a_gamma;
//                 Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _pip_mu_prime_cm_.Vect().Unit();
//
//                 a_beta = sqrt(1. / (1 - pow((_pip_mu_prime_cm_.Vect().Unit() * _pim_mu_prime_cm_.Vect().Unit()),
//                 2))); b_beta = -(_pip_mu_prime_cm_.Vect().Unit() * _pim_mu_prime_cm_.Vect().Unit()) * a_beta;
//                 Vect3_beta = a_beta * _pim_mu_prime_cm_.Vect().Unit() + b_beta * _pip_mu_prime_cm_.Vect().Unit();
//
//                 alpha_PPIm_piPIp = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
//
//                 if (Vect3_gamma.Cross(Vect3_beta) * _pip_mu_prime_cm_.Vect() < 0)
//                         alpha_PPIm_piPIp = 360. - alpha_PPIm_piPIp; // α[pπ−][p'π+]
//                 //  return;
//                 // return (phi_P_cm, phi_PIp_cm, phi_PIm_cm, theta_PIm_cm,theta_PIp_cm,
//                 // theta_P_cm, alpha_PPIp_piPIm,alpha_PIpPIm_pipf,alpha_PPIm_piPIp);
//
//                 _theta_gamma = theta_gamma;
//                 _phi_gamma = phi_gamma;
//                 _theta_prot = theta_prot;
//                 _phi_prot = phi_prot;
//                 _theta_pip = theta_pip;
//                 _phi_pip = phi_pip;
//                 _theta_pim = theta_pim;
//                 _phi_pim = phi_pim;
//
//                 _alpha_ppip_pipim = alpha_PPIp_piPIm;
//                 _alpha_pippim_pipf = alpha_PIpPIm_pipf; // i haven't made functions to output these angles
//                 _alpha_ppim_pipip = alpha_PPIm_piPIp; // i haven't made functions to output these angles
//         }
// }
// float Reaction::alpha_ppip_pipim() {
//         return _alpha_ppip_pipim;
// }

float Reaction::MM() { return _MM; }
float Reaction::MM2() { return _MM2; }
float Reaction::MM_wop() { return _MM_wop; }
float Reaction::MM2_wop() { return _MM2_wop; }

float Reaction::W() { return _W; }
float Reaction::Q2() { return _Q2; }

float Reaction::beta() { return _beta; }
float Reaction::gamma_() { return _gamma; }

float Reaction::W_ep() { return _W_ep; }
float Reaction::W_2pi() { return _W_2pi; }
float Reaction::W_delta_pp() { return _W_delta_pp; }
float Reaction::W_delta_zero() { return _W_delta_zero; }
float Reaction::W_rho() { return _W_rho; }

float Reaction::W_singlepip() { return _W_singlepip; }

float Reaction::Q2_2pi() { return _Q2_2pi; }
