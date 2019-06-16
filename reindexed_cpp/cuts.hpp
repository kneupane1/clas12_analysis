
#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD

#include "constants.hpp"
class Cuts {
 private:
  bool _good_e;
  bool _good_p;
  bool _good_pip;
  bool _good_pim;
  bool _good_hadron;

  float x_PCAL_rot, y_PCAL_rot, angle, height_PCAL, slope_PCAL, left_PCAL, right_PCAL, radius2_PCAL, dcR1_height,
      dcR2_height, dcR3_height, x1_rot, y1_rot, x2_rot, y2_rot, x3_rot, y3_rot, slope, left_r1, right_r1, left_r2,
      right_r2, left_r3, right_r3, radius2_DCr1, radius2_DCr2, radius2_DCr3;

 public:
  Cuts();
  ~Cuts();
  bool electron_cuts(int status, int charge, float sf, float vertex_pos, float chi_sq, float mom_el, float th_el,
                     float ph_el, int sec, float x_PCAL, float y_PCAL, float dc_r1x, float dc_r1y, float dc_r2x,
                     float dc_r2y, float dc_r3x, float dc_r3y);
  bool proton_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
  bool hadron_cuts_ctof(int status, int charge, float mom, int pid, float chi_sq);
  bool pip_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
  bool pim_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
};

#endif
