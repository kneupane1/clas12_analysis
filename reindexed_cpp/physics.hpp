/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD
#include "TROOT.h"
#include "constants.hpp"
#include <TLorentzVector.h>
#include <TRotation.h>
namespace physics {
// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime);
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime);
double xb_calc(double Q2, double E_prime);
// overload with 4 vectors instaed of other calculations
double xb_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime);
double theta_calc(double cosz);
double phi_calc(double cosx, double cosy);

// boost to COM frame
float q_3(TLorentzVector e_mu, TLorentzVector e_mu_prime);
float_t beta_boost(TLorentzVector e_mu, TLorentzVector e_mu_prime);
double gamma(TLorentzVector e_mu, TLorentzVector e_mu_prime);
TLorentzVector boost_(TLorentzVector four_vect, TLorentzVector e_mu,
                      TLorentzVector e_mu_prime);

double theta_fn(TLorentzVector four_vect, TLorentzVector e_mu,
                TLorentzVector e_mu_prime);

// void TLorentzVector Boost(Double_t bx, Double_t by, Double_t bz);
// namespace physics
// double gamma(TLorentzVector e_mu) {
//   return ((e_mu.E() + MASS_P * pow(c_special_units, 2)) /
//           (((MASS_E * MASS_E + MASS_P * MASS_P) * pow(c_special_units, 4)) +
//            2 * e_mu.E() * MASS_P * pow(c_special_units, 2)));
// }
//
// TVector3 beta_boost(TLorentzVector e_mu) {
//   return (e_mu_3 * c_special_units /
//           (e_mu.E() + MASS_P * pow(c_special_units, 2)));
// }
//
// TLorentzVector P_eCM(TLorentzVector e_mu) {
//   // double Energy_eCM = gamma* ( e_mu.E() -
//   // (beta_boost+e_mu_3).Mag()*c_special_units;
//   TVector3 P_e_CM_3 =
//       gamma * (e_mu_3 - beta_boost * e_mu.E() * c_special_units);
//   TLorentzVector e_mu_cm;
//   e_mu_cm.SetVectM(P_e_CM_3, Energy_eCM);
//   return e_mu_cm;
}
#endif
