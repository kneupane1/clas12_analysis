
#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "TMath.h"
#include "constants.hpp"
#include "histogram.hpp"
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
                     float ph_el, int sec, float x_PCAL, float y_PCAL, float dc_r1x, float dc_r1y/*, float dc_r2x,
                     float dc_r2y, float dc_r3x, float dc_r3y*/);
  bool proton_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
  bool pip_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
  bool pim_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
};

#endif
