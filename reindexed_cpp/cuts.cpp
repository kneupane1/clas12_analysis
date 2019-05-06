#include <iostream>
#include "cuts.hpp"
Cuts::Cuts() {
  _status = std::nan("-99");
  _charge = std::nan("-99");
  _min_mom = std::nan("-99");
  _sf = std::nan("-99");
  _vertex_pos = std::nan("-99");

  _good_e = false;
}
Cuts::~Cuts() {}

bool Cuts::electron_cuts(int status, int charge, float min_mom, float sf, float vertex_) {
  if (2000 < status < 4000) {
    if (charge == -1) {
      if (min_mom > 1.0) {
        if (sf > 0.15 && sf < 0.28) {
          if (-10.0 < vertex_p < 5.0) {
            _good_e = true;
          }
        }
      }
    }
  }
  return _good_e;
}
// bool Cuts::ELE(){
//         return _good_e;
// }
