#include "electron_cuts.hpp"
cuts::cuts() {
  _charge = std::nan("-99");
  _dc = std::nan("-99");
  _cc = std::nan("-99");
  _ec = std::nan("-99");
  _sc = std::nan("-99");
  _dc_sector = std::nan("-99");
  _good_e = kFALSE;
}
cuts::~cuts() {}

bool electron_cuts(int charge, int dc, int cc, int sc, int ec, float sf) {
  if (charge == -1) {
    if (dc > 0) {
      if (cc > 0) {
        if (ec > 0) {
          if (sc > 0) {
            if (sf > 0.15 && sf < 0.28) {
              _good_e = kTRUE;
            }
          }
        }
      }
    }
  }
  return _good_e;
}
