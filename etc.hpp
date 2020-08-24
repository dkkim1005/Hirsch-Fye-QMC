#pragma once

#include <fftw3.h>
#include <complex>
#include <vector>
#include <cmath>
#include <string>

namespace etc
{
// FFTW3 library.
// FFTW_DIRECTION : FORWARD(-1) , BACKWARD(+1)
enum FOURIER_TRANSFORM { FORWARD = -1, BACKWARD = 1 };

// Fast Fourier transform for 1-D.
inline void fft_1d_complex(const std::complex<double> * input, std::complex<double> * output,
  const int N, const FOURIER_TRANSFORM FFTW_DIRECTION)
{
  fftw_complex * in = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*N)),
    * out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*N));
  fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_DIRECTION, FFTW_ESTIMATE);

  for(int i=0; i<N; ++i)
  {
    in[i][0] = input[i].real();
    in[i][1] = input[i].imag();
  }   

  fftw_execute(p); // repeat as needed 

  for(int i=0; i<N; ++i)
    output[i] = std::complex<double>(out[i][0], out[i][1]);

  fftw_destroy_plan(p);

  fftw_free(in);
  fftw_free(out);
}

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

inline std::string to_string(const double & value)
{
  std::string valstr = std::to_string(value);
  valstr.erase(valstr.find_last_not_of('0') + 1, std::string::npos);
  valstr.erase(valstr.find_last_not_of('.') + 1, std::string::npos);
  return valstr;
}
} // end namespace etc
