#pragma once

#include <vector>
#include <cmath>
#include <complex>
#include <functional>
#include "etc.hpp"

typedef std::vector<double> ImTimeGreen;
typedef std::vector<std::complex<double> > ImFreqGreen;

// Inverse Fourier transform of the fermionic Matsubara green function (frequency -> time).
void GreenInverseFourier(const ImFreqGreen & giwn, const double beta,
  const std::vector<double> & M, ImTimeGreen & gtau);

// Fourier transform of the fermionic Matsubara green function (time -> frequency).
void GreenFourier(const ImTimeGreen & gtau, const double beta, ImFreqGreen & giwn, const bool useSpline = true);

// D : half bandwidth, ne : mesh # of energy domain
ImFreqGreen SemiCircularGreen(const int niwn, const double beta,
  const double mu, const double V, const double D, const int ne = 1001);
