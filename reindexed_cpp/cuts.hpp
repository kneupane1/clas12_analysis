
#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "constants.hpp"
class Cuts {
 private:
  int _status;
  int _charge;
  float _min_mom;
  float _sf;
  float _vertex_pos;
  bool _good_e;

 public:
  Cuts();
  ~Cuts();
  bool electron_cuts(int status, int charge, float min_mom, float sf, float vertex_);
  bool hadron_cuts();
};

#endif
