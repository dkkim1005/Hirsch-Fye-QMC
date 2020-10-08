#pragma once

// Ref: 1) https://www.physics.rutgers.edu/grad/509/qmc.pdf
//      2) Hirsch-Fye QMC solver, ALPS project (http://alps.comp-phys.org/mediawiki/index.php/Main_Page)

//#define USE_TRNG4

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#ifdef USE_TRNG4 
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#else
#include <random>
#endif
#include "blas_lapack.hpp"
#include "etc.hpp"

typedef std::vector<double> ImTimeGreen;

#ifdef USE_TRNG4
template <typename T>
using uniform_dist = trng::uniform01_dist<T>;
using rand_gen = trng::yarn2;
#else
template <typename T>
using uniform_dist = std::uniform_real_distribution<T>;
using rand_gen = std::mt19937;
#endif

class HirschFyeQMC
{
public:
  HirschFyeQMC(const double * gtau0, const int ntau, const double beta, const double U, const unsigned long seed);
  // return acceptance rate
  double do_montecarlo_step(const int nmcsteps, const int nperiodSweeps);
  void accumulate(double * G_up, double * G_dw, double * G2_up, double * G2_dw);
  void get_auxiliary_fields(int * spins) const
  {
    std::memcpy(spins, spins_.data(), sizeof(int)*spins_.size());
  }

private:
  enum SPIN_STATE {UP = 1, DOWN = -1};
  double e1lamb_(const int & s) const { return exp1lamb_[s+1]; }
  double e2lamb_(const int & s) const { return exp2lamb_[s+1]; }
  void clean_update_(std::vector<double> & gtaumat, SPIN_STATE SPIN);
  bool single_flip_update_();
  void green_vector_from_matrix_(std::vector<double> & gtauvec, const std::vector<double> & gtaumat);

  int nsweeps_, nchain_;
  double acceptanceRate_;
  const int N_;
  const double lambda_;
  const std::vector<double> gtau0_;
  std::vector<double> gtau0mat_, gtau_up_, gtau_dw_, Gtau_up_, Gtau_dw_, vi_up_, uj_up_, vi_dw_, uj_dw_, e_ls_, A_;
  std::vector<int> spins_;
  uniform_dist<double> unidist_;
  rand_gen randdev_;
  std::array<double, 3> exp1lamb_, exp2lamb_;
};
