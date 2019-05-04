
#ifndef ELECTRON_CUTS_H_GUARD
#define ELECTRON_CUTS_H_GUARD
#include <iostream>
#include "constants.hpp"

class cuts {
 private:
  int _charge, _dc, _sc, _cc, _ec;
  float sf;
  int dc_sectors;
  int _status;
  bool _good_e;

 public:
  electron_cuts();
  ~electron_cuts();

  bool electron_cuts(int charge, int status, int dc, int sc, int cc, int ec, float sf, int dc_sectors);
  bool hadron_cuts();
};

#endif
