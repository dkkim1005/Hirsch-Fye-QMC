#pragma once

#include <fftw3.h>
#include <complex>
#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <assert.h>

namespace etc
{
// FFTW3 library.
// FFTW_DIRECTION : FORWARD(-1) , BACKWARD(+1)
enum FOURIER_TRANSFORM { FORWARD = -1, BACKWARD = 1 };

// Fast Fourier transform for 1-D.
void fft_1d_complex(const std::complex<double> * input, std::complex<double> * output,
  const int N, const FOURIER_TRANSFORM FFTW_DIRECTION);

// Ref: https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_rule
template <typename FuncType>
typename FuncType::returnType simpson_rule(std::vector<double> & x, const FuncType & func)
{
  const int N = x.size()-1;
  std::vector<double> h(N, 0);
  std::vector<typename FuncType::returnType> y(x.size(), 0);
  for (int i=0; i<N; ++i)
    h[i] = x[i+1] - x[i];
  for (int i=0; i<x.size(); ++i)
    y[i] = func(x[i]);
  typename FuncType::returnType result = 0;
  for (int i=1; i<N; i+=2)
  {
    const double hph = h[i] + h[i-1];
    result += y[i]*(std::pow(h[i], 3) + std::pow(h[i-1], 3) + 3.0*h[i]*h[i-1]*hph)/(6.0*h[i]*h[i-1]);
    result += y[i-1]*(2.0*std::pow(h[i-1], 3) - std::pow(h[i], 3) + 3.0*h[i]*std::pow(h[i-1], 2))/(6.0*h[i-1]*hph);
    result += y[i+1]*(2.0*std::pow(h[i], 3) - std::pow(h[i-1], 3) + 3.0*h[i-1]*std::pow(h[i], 2))/(6.0*h[i]*hph);
  }
  if ((N+1)%2 == 0)
  {
    result += y[N]*(2*std::pow(h[N-1], 2) + 3.0*h[N-2]*h[N-1])/(6.0*(h[N-2]+h[N-1]));
    result += y[N-1]*(std::pow(h[N-1], 2) + 3.0*h[N-1]*h[N-2])/(6.0*h[N-2]);
    result -= y[N-2]*std::pow(h[N-1], 3)/(6.0*h[N-2]*(h[N-2]+h[N-1]));
  }
  return result;
}

template <typename T>
void transpose(T * mat, const int N)
{
  for (int i=0; i<N; ++i)
    for (int j=i+1; j<N; ++j)
    {
      const T temp = mat[i*N+j];
      mat[i*N+j] = mat[j*N+i];
      mat[j*N+i] = temp;
    }
}

std::string to_string(const double & value);

class CubicHermiteSpline
{
public:
  CubicHermiteSpline(const double* xVec, const double* yVec, const double* mVec, const size_t size);
  virtual ~CubicHermiteSpline() { this->delete_mem_(); }
  double operator()(const double x) const;

protected:
  double * xVec_, * yVec_, * mVec_;
  size_t size_;

  explicit CubicHermiteSpline();
  void delete_mem_();
  void allocate_mem_(const double* xVec, const double* yVec, const double* mVec, const size_t size);
  size_t binary_search_(const double& x) const;
  double interp_f_(const double& t, const size_t& idx) const;
};

class MonotoneCubicInterpolation : public CubicHermiteSpline
{
public:
  MonotoneCubicInterpolation(const double* xVec, const double* yVec, const size_t size);
  virtual ~MonotoneCubicInterpolation() {}
};
} // end namespace etc
