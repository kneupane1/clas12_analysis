
#ifndef ELECTRON_CUTS_H_GUARD
#define ELECTRON_CUTS_H_GUARD
#include <iostream>
#include "constants.hpp"

class cuts {
 private:
  int _status;
  int _charge;
  float _min_mom;
  float sf;
  float _vertex_pos;
  bool _good_e;

 public:
  electron_cuts();
  ~electron_cuts();

  bool electron_cuts(int status, int charge, float min_mom, float sf, float vertex_);
  bool hadron_cuts();
};

#endif
