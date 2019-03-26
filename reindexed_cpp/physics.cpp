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

float_t beta_boost(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q = e_mu - e_mu_prime;
  return ((sqrt(q[3] * q[3] + Q2_calc(e_mu, e_mu_prime))) / (q[3] + MASS_P));
}
double gamma(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  double b = beta_boost(e_mu, e_mu_prime);
  return sqrt(1 / (b * b));
}
float q_3(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q = e_mu - e_mu_prime;
  return Q2_calc(e_mu, e_mu_prime);
}

TLorentzVector boost_(TLorentzVector four_vect, TLorentzVector e_mu,
                      TLorentzVector e_mu_prime) {
  TRotation rot;
  TLorentzVector q = e_mu - e_mu_prime;
  //  double beta = sqrt((q[3] * q[3] + Q2_calc(e_mu, e_mu_prime)) /
  //                   (q[3] + MASS_P)); // q[3] = q.Vect().Mag()
  float_t beta_1 = beta_boost(e_mu, e_mu_prime);
  TVector3 uz = q.Vect().Unit(); // uit vector along virtual photon
  TVector3 ux = (e_mu.Vect().Cross(e_mu_prime.Vect()))
                    .Unit(); // unit vector along e cross e'
  ux.Rotate(3. * PI / 2,
            uz); // rotating ux by 3pi/2 with uz as axis of roration
  rot.SetZAxis(uz, ux).Invert(); // setting TRotation rot
  four_vect.Transform(rot);
  four_vect.Boost(0, 0, -beta_1); // -beta ko value (0.5 to -0.5 huda samma
                                  // value aauchha nattra aaudyna)
  // four_vect.Boost(-(e_mu + target).BoostVector());
  return four_vect;
}

double theta_fn(TLorentzVector four_vect, TLorentzVector e_mu,
                TLorentzVector e_mu_prime) {
  TLorentzVector q = e_mu - e_mu_prime;
  return acos((q.Vect().Dot(four_vect.Vect())) /
              (q.Vect().Mag() * (four_vect.Vect()).Mag()));
}
// TLorentzVector P_eCM(TLorentzVector e_mu) {
//   // double Energy_eCM = gamma* ( e_mu.E() -
//   // (beta_boost+e_mu_3).Mag()*c_special_units;
//   TVector3 P_e_CM_3 =
//       gamma * (e_mu_3 - beta_boost * e_mu.E() * c_special_units);
//   TLorentzVector e_mu_cm;
//   e_mu_cm.SetVectM(P_e_CM_3, Energy_eCM);
//   return e_mu_cm;
// }
}
// namespace physics
