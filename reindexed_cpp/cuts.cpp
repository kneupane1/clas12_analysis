
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
bool Cuts::electron_cuts(int status, int charge, float sf, float vertex_pos, float chi_sq, float mom_el, float th_el,
                         float ph_el) {
  if (2000 <= status && status < 4000) {
    if (charge == -1) {
      if (mom_el > 1.0) {
        if (sf > 0.18 && sf < 0.28) {
          if (-10.0 < vertex_pos && vertex_pos < 5.0) {
            if (-2000 < chi_sq && chi_sq < 2000) {
              th_min = (11.7398 + 8.21504 / (0.433327 * mom_el + 0.158076));
              par1 = 0.85 + 1.1 * mom_el;
              par2 = -62.8 - 30. * mom_el;
              par3 = 0.0047 * mom_el + 0.0079;

              fid_a = 41.3 * pow((sin((th_el - th_min) * par3)), (par1 + par2 / th_el + 1485. / th_el / th_el)) + 1.;
              fid_b = -41.3 * pow((sin((th_el - th_min) * par3)), (par1 + par2 / th_el + 1485. / th_el / th_el)) - 1.;

              th_max = 76.8617 - 76.537 * mom_el + 77.9387 * mom_el * mom_el - 28.389 * mom_el * mom_el * mom_el;

              if ((th_el > th_min) && (th_el < th_max) && (ph_el > fid_b + 360) && (ph_el < fid_a + 360)) {
                //  if ((mom_el < 1.75999) && (mom_el > 0.4)) {
                _good_e = true;
                //  }
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
