/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DT_H_GUARD
#define DT_H_GUARD
#include "constants.hpp"
#include <iostream>

class Delta_T {
private:
  std::vector<double> masses = {MASS_E, MASS_P, MASS_PIP,
                                MASS_KP}; // std::vector is a sequence
                                          // container that encapsulates dynamic
                                          // size arrays.
  float _sc_t_v = std::nanf("-99");
  float _sc_r_v = std::nanf("-99");
  float _vertex = std::nanf("-99");
  float _sc_t = std::nanf("-99");
  float _sc_r = std::nanf("-99");
  float _sc_t_ctof = std::nanf("-99"); /// maile thapeko
  float _sc_r_ctof = std::nanf("-99"); // maile thapeko

  float _beta = std::nanf("-99");
  float _momentum = std::nanf("-99");

  float _vertex_time(float sc_time, float sc_pathlength,
                     float relatavistic_beta);
  float _deltat(int num);
  float _dt_calc_ctof(int num);

public:
  Delta_T(float time_1b, float path_1b, float time_1a, float path_1a,
          float time_2, float path_2);
  ~Delta_T();

  void dt_calc(float momentum, float time_1b, float path_1b, float time_1a,
               float path_1a, float time_2, float path_2, float time_ctof,
               float path_ctof);
  float dt_E();
  float dt_P();
  float dt_Pi();
  float dt_K();
  float momentum();

  void dt_calc_1(float momentum, float time_ctof, float path_ctof);

  float dt_ctof_E();
  float dt_ctof_P();
  float dt_ctof_Pi();
  float dt_ctof_K();

  /*
  TODO:
  Add Parts for each detector 1a/1b/2/ctof

  float dt_1b_E(int part);
  float dt_1b_P(int part);
  float dt_1b_Pi(int part);
  float dt_1b_K(int part);

  float dt_1a_E(int part);
  float dt_1a_P(int part);
  float dt_1a_Pi(int part);
  float dt_1a_K(int part);
  etc...
  */
};

#endif
