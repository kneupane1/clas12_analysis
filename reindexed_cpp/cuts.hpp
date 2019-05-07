
#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "constants.hpp"
class Cuts {
private:

bool _good_e;
bool _good_p;
bool _good_pip;
bool _good_pim;

public:
Cuts();
~Cuts();
bool electron_cuts(int status, int charge, float min_mom, float sf, float vertex_pos,float chi_sq);
bool proton_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
bool pip_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
bool pim_cuts(int status, int charge, float min_mom, int pid, float chi_sq);
};

#endif
