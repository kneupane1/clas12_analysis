
#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "TMath"
#include "constants.hpp"
#include "reaction.hpp"

class Cuts {
 private:
  bool _good_e;
  bool _good_p;
  bool _good_pip;
  bool _good_pim;
  Float_t th_min, th_max, par1, par2, par3, fid_a, fid_b;

 public:
  Cuts();
  ~Cuts();
  bool electron_cuts(int status, int charge, float sf, float vertex_pos, float chi_sq, float mom_el, float th_el,
                     float ph_el);
  bool proton_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
  bool pip_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
  bool pim_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
};

#endif
