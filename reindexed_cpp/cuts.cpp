
#include "cuts.hpp"
using namespace std;
Cuts::Cuts() {
  _status = std::nan("-99");
  _charge = std::nan("-99");
  _min_mom = std::nan("-99");
  _sf = std::nan("-99");
  _vertex_pos = std::nan("-99");

  _good_e = false;
}
Cuts::~Cuts() {}
bool Cuts::electron_cuts(int status, int charge, float min_mom, float sf, float vertex_pos) {
  if (2000 < status && status < 4000) {
    if (charge == -1) {
      if (min_mom > 1.0) {
        if (sf > 0.18 && sf < 0.28) {
          if (-10.0 < vertex_pos && vertex_pos < 5.0) {
            _good_e = true;
          }
        }
      }
    }
  }
  return _good_e;
}
