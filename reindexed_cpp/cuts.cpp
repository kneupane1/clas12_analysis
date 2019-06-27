
#include "cuts.hpp"

Cuts::Cuts() {
  _good_e = false;
  _good_p = false;
  _good_pip = false;
  _good_pim = false;
  _good_hadron = false;

  float x_PCAL_rot, y_PCAL_rot, angle, height_PCAL, slope_PCAL, left_PCAL, right_PCAL, radius2_PCAL, dcR1_height,
      dcR2_height, dcR3_height, x1_rot, y1_rot, x2_rot, y2_rot, x3_rot, y3_rot, slope, left_r1, right_r1, left_r2,
      right_r2, left_r3, right_r3, radius2_DCr1, radius2_DCr2, radius2_DCr3;
}
Cuts::~Cuts() {}

bool Cuts::electron_cuts(int status, int charge, float sf, float vertex_pos, float chi_sq, float mom_el, float th_el,
                         float ph_el, int sec, float x_PCAL, float y_PCAL, float dc_r1x, float dc_r1y, float dc_r2x,
                         float dc_r2y, float dc_r3x, float dc_r3y) {
  if (2000 <= status && status < 4000) {
    if (charge == -1) {
      if (mom_el > 0.20) {
        if (sf > 0.18 && sf < 0.28) {
          if (-10.0 < vertex_pos && vertex_pos < 5.0) {
            if (-2000 < chi_sq && chi_sq < 2000) {
              x_PCAL_rot = y_PCAL * sin(sec * 60.0 * PI / 180) + x_PCAL * cos(sec * 60.0 * PI / 180);
              y_PCAL_rot = y_PCAL * cos(sec * 60.0 * PI / 180) - x_PCAL * sin(sec * 60.0 * PI / 180);
              angle = 60;
              height_PCAL = 45;
              slope_PCAL = 1 / tan(0.5 * angle * PI / 180);
              left_PCAL = (height_PCAL - slope_PCAL * y_PCAL_rot);
              right_PCAL = (height_PCAL + slope_PCAL * y_PCAL_rot);
              radius2_PCAL = pow(height_PCAL + 6, 2) - pow(y_PCAL_rot, 2);

              if (x_PCAL_rot > left_PCAL && x_PCAL_rot > right_PCAL && pow(x_PCAL_rot, 2) > radius2_PCAL &&
                  x_PCAL_rot < 372) {
                angle = 60;
                dcR1_height = 15;  // 31
                dcR2_height = 25;  // 47
                dcR3_height = 35;  // 53

                x1_rot = dc_r1y * sin(sec * 60.0 * PI / 180) + dc_r1x * cos(sec * 60.0 * PI / 180);
                y1_rot = dc_r1y * cos(sec * 60.0 * PI / 180) - dc_r1x * sin(sec * 60.0 * PI / 180);
                x2_rot = dc_r2y * sin(sec * 60.0 * PI / 180) + dc_r2x * cos(sec * 60.0 * PI / 180);
                y2_rot = dc_r2y * cos(sec * 60.0 * PI / 180) - dc_r2x * sin(sec * 60.0 * PI / 180);
                x3_rot = dc_r3y * sin(sec * 60.0 * PI / 180) + dc_r3x * cos(sec * 60.0 * PI / 180);
                y3_rot = dc_r3y * cos(sec * 60.0 * PI / 180) - dc_r3x * sin(sec * 60.0 * PI / 180);

                slope = 1 / tan(0.5 * angle * PI / 180);

                left_r1 = (dcR1_height - slope * y1_rot);
                right_r1 = (dcR1_height + slope * y1_rot);
                left_r2 = (dcR2_height - slope * y2_rot);
                right_r2 = (dcR2_height + slope * y2_rot);
                left_r3 = (dcR3_height - slope * y3_rot);
                right_r3 = (dcR3_height + slope * y3_rot);

                radius2_DCr1 = pow(15, 2) - pow(y1_rot, 2);  // 32 stafen
                radius2_DCr2 = pow(20, 2) - pow(y2_rot, 2);  // 49
                radius2_DCr3 = pow(25, 2) - pow(y3_rot, 2);  // 54

                if (x1_rot > left_r1 && x1_rot > right_r1 && pow(x1_rot, 2) > radius2_DCr1) {
                  if (x2_rot > left_r2 && x2_rot > right_r2 && pow(x2_rot, 2) > radius2_DCr2) {
                    if (x3_rot > left_r3 && x3_rot > right_r3 && pow(x3_rot, 2) > radius2_DCr3) {
                      _good_e = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return _good_e;
}
bool Cuts::hadron_cuts_ctof(int status, int charge, float mom, int pid, float chi_sq) {
  if (4000 <= status && status <= 6000) {
    if (charge != 0) {
      if (mom > 0.20 && mom < 3.0) {
        // if (pid == 2212 || pid == 211 || pid == -211) {
        if (-2000 < chi_sq && chi_sq < 2000) {
          _good_hadron = true;
        }
        //}
      }
    }
  }
  return _good_hadron;
}
bool Cuts::proton_cuts(int status, int charge, float min_mom, int pid, float chi_sq) {
  if (2000 <= status && status <= 4000) {  // forward ko lagi 2000 to 4000 and central ko lagi >= 4000
    if (charge != 0) {
      if (min_mom > 0.30) {
        if (pid == 2212) {
          if (-2000 < chi_sq && chi_sq < 2000) {
            _good_p = true;
          }
        }
      }
    }
  }
  return _good_p;
}

bool Cuts::pip_cuts(int status, int charge, float min_mom, int pid, float chi_sq) {
  if (2000 <= status && status <= 4000) {
    if (charge != 0) {
      if (min_mom > 0.30) {
        if (pid == 211) {
          if (-2000 < chi_sq && chi_sq < 2000) {
            _good_pip = true;
          }
        }
      }
    }
  }
  return _good_pip;
}
bool Cuts::pim_cuts(int status, int charge, float min_mom, int pid, float chi_sq) {
  if (2000 <= status && status <= 4000) {
    if (charge != 0) {
      if (min_mom > 0.30) {
        if (pid == -211) {
          if (-2000 < chi_sq && chi_sq < 2000) {
            _good_pim = true;
          }
        }
      }
    }
  }
  return _good_pim;
}
