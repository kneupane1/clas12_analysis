
#include "cuts.hpp"
using namespace std;
Cuts::Cuts() {
  _good_e = false;
  _good_p = false;
  _good_pip = false;
  _good_pim = false;
  th_min = std::nan("-99");
  th_max = std::nan("-99");
  par1 = std::nan("-99");
  par2 = std::nan("-99");
  par3 = std::nan("-99");
  fid_a = std::nan("-99");
  fid_b = std::nan("-99");
}
Cuts::~Cuts() {}
Histogram *hist1 = new Histogram();

bool Cuts::electron_cuts(int status, int charge, float sf, float vertex_pos, float chi_sq, float mom_el, float th_el,
                         float ph_el, int sec, float x_PCAL, float y_PCAL, float dc_r1x, float dc_r1y, float dc_r2x,
                         float dc_r2y, float dc_r3x, float dc_r3y) {
  if (2000 <= status && status < 4000) {
    if (charge == -1) {
      if (mom_el > 0.30) {
        if (sf > 0.18 && sf < 0.28) {
          if (-10.0 < vertex_pos && vertex_pos < 5.0) {
            if (-2000 < chi_sq && chi_sq < 2000) {
              float x_PCAL_rot = y_PCAL * sin(sec * 60.0 * PI / 180) + x_PCAL * cos(sec * 60.0 * PI / 180);
              float y_PCAL_rot = y_PCAL * cos(sec * 60.0 * PI / 180) - x_PCAL * sin(sec * 60.0 * PI / 180);
              float angle_PCAL = 60;
              float height_PCAL = 45;
              float slope_PCAL = 1 / tan(0.5 * angle_PCAL * PI / 180);
              float left_PCAL = (height_PCAL - slope_PCAL * y_PCAL_rot);
              float right_PCAL = (height_PCAL + slope_PCAL * y_PCAL_rot);
              float radius2_PCAL = pow(height_PCAL + 6, 2) - pow(y_PCAL_rot, 2);

              if (x_PCAL_rot > left_PCAL && x_PCAL_rot > right_PCAL && pow(x_PCAL_rot, 2) > radius2_PCAL &&
                  x_PCAL_rot < 372) {
                float dc_angle = 60;
                float dcR1_height = 31;
                float dcR2_height = 47;
                float dcR3_height = 53;

                float x1_rot = dc_r1x * sin(sec * 60.0 * Pival / 180) + dc_r1x * cos(sec * 60.0 * Pival / 180);
                float y1_rot = dc_r1y * cos(sec * 60.0 * Pival / 180) - dc_r1y * sin(sec * 60.0 * Pival / 180);
                float x2_rot = dc_r2x * sin(sec * 60.0 * Pival / 180) + dc_r2x * cos(sec * 60.0 * Pival / 180);
                float y2_rot = dc_r2y * cos(sec * 60.0 * Pival / 180) - dc_r2y * sin(sec * 60.0 * Pival / 180);
                float x3_rot = dc_r3x * sin(sec * 60.0 * Pival / 180) + dc_r3x * cos(sec * 60.0 * Pival / 180);
                float y3_rot = dc_r3y * cos(sec * 60.0 * Pival / 180) - dc_r3y * sin(sec * 60.0 * Pival / 180);

                float slope = 1 / tan(0.5 * dc_angle * Pival / 180);

                float left_r1 = (height - slope * y1_rot);
                float right_r1 = (height + slope * y1_rot);
                float left_r2 = (height - slope * y2_rot);
                float right_r2 = (height + slope * y2_rot);
                float left_r3 = (height - slope * y3_rot);
                float right_r3 = (height + slope * y3_rot);

                float radius2_DCr1 = pow(32, 2) - pow(y1_rot, 2);
                float radius2_DCr2 = pow(49, 2) - pow(y2_rot, 2);
                float radius2_DCr3 = pow(54, 2) - pow(y3_rot, 2);

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
bool Cuts::proton_cuts(int status, int charge, float min_mom, int pid, float chi_sq) {
  if (2000 <= status && status < 6000) {  // forward ko lagi 2000 to 4000 and central ko lagi >= 4000
    if (charge != 0) {
      if (min_mom > 0.20) {
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
  if (2000 <= status && status < 6000) {
    if (charge != 0) {
      if (min_mom > 0.20) {
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
  if (2000 <= status && status < 6000) {
    if (charge != 0) {
      if (min_mom > 0.20) {
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
