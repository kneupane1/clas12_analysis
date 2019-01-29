
#ifndef ELECTRON_CUTS_H_GUARD
#define ELECTRON_CUTS_H_GUARD
#include "constants.hpp"
#include <iostream>

class cuts {
private:
  int _charge, _dc, _sc, _cc, _ec;
  float sf;
  int dc_sectors;
  bool _good_e;

public:
  electron_cuts();
  ~electron_cuts();

  bool electron_cuts(int charge, int dc, int sc, int cc, int ec, float sf,
                     int dc_sectors);
  bool hadron_cuts();
};

#endif
