/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
#include "physics.hpp"

namespace physics {
// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  return -q_mu.Mag2();
}
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  TVector3 p_mu_3(0, 0, 0);
  TLorentzVector p_mu;
  p_mu.SetVectM(p_mu_3, MASS_P);
  return (p_mu + q_mu).Mag();
}

double xb_calc(double Q2, double E_prime) {
  double gamma = CLAS12_E - E_prime;
  double xb = (Q2 / (2 * MASS_P * gamma));
  return xb;
}
// overload with 4 vectors
double xb_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  double Q2 = Q2_calc(e_mu, e_mu_prime);
  TLorentzVector q = e_mu - e_mu_prime;
  TLorentzVector target(0, 0, 0, MASS_P);
  return (Q2 / (2 * (q.Dot(target))));
}
// double vertex_time_ctof(double sc_time, double sc_pathlength,
//                         double relatavistic_beta) {
//   return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
// }
//
// double deltat_ctof(double electron_vertex_time, double mass, double momentum,
//                    double sc_t, double sc_r) {
//   double relatavistic_beta =
//       1.0 / sqrt(1.0 + (mass / momentum) * (mass / momentum));
//   return electron_vertex_time - vertex_time_ctof(sc_t, sc_r,
//   relatavistic_beta);
// }
}
// namespace physics
