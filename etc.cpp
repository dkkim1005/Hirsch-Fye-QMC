#include "etc.hpp"

namespace etc
{
void fft_1d_complex(const std::complex<double> * input, std::complex<double> * output,
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

std::string to_string(const double & value)
{
  std::string valstr = std::to_string(value);
  valstr.erase(valstr.find_last_not_of('0') + 1, std::string::npos);
  valstr.erase(valstr.find_last_not_of('.') + 1, std::string::npos);
  return valstr;
}


CubicHermiteSpline::CubicHermiteSpline(const double* xVec, const double* yVec, const double* mVec, const size_t size):
  xVec_(nullptr),
  yVec_(nullptr),
  mVec_(nullptr),
  size_(0) { this->allocate_mem_(xVec, yVec, mVec, size); }

double CubicHermiteSpline::operator()(const double x) const
{
  const size_t idx = this->binary_search_(x);
  if(idx == size_ - 1)
    return yVec_[size_-1];
  const double t = (x - xVec_[idx])/(xVec_[idx+1] - xVec_[idx]);
  return this->interp_f_(t, idx);
}

CubicHermiteSpline::CubicHermiteSpline():
  xVec_(nullptr),
  yVec_(nullptr),
  mVec_(nullptr),
  size_(0) {}

void CubicHermiteSpline::delete_mem_()
{
  if(xVec_ != nullptr)
    delete [] xVec_;
  if(yVec_ != nullptr)
    delete [] yVec_;
  if(mVec_ != nullptr)
    delete [] mVec_;
}

void CubicHermiteSpline::allocate_mem_(const double* xVec, const double* yVec, const double* mVec, const size_t size)
{
  this->delete_mem_();
  xVec_ = new double [size];
  yVec_ = new double [size];
  mVec_ = new double [size];
  size_ = size;
  std::memcpy(xVec_, xVec,sizeof(double)*size);
  std::memcpy(yVec_, yVec,sizeof(double)*size);
  std::memcpy(mVec_, mVec,sizeof(double)*size);
}

size_t CubicHermiteSpline::binary_search_(const double& x) const
{
  assert(xVec_[0] <= x and x <= xVec_[size_-1]);
  size_t idx_l = 0, idx_r = size_ - 1, idx = size_/2;
  auto is_x_in_Boundary = [this, &x](const size_t& idx) -> bool
    { return this -> xVec_[idx] <= x and x < this -> xVec_[idx + 1];};
  while (1)
  {
    if (idx_r - idx_l == 1) 
    {
      if (is_x_in_Boundary(idx))
        return idx;
      else
        return idx + 1;
    }
    if (is_x_in_Boundary(idx))
      return idx;
    else if (xVec_[idx+1] <= x)
    {
      idx_l = idx;
      idx = (idx_r - idx_l)/2 + idx_l;
    }
    else
    {
      idx_r = idx;
      idx = (idx_r - idx_l)/2 + idx_l;
    }
  }
}

double CubicHermiteSpline::interp_f_(const double& t, const size_t& idx) const
{
  return (2*std::pow(t,3) - 3*std::pow(t,2) + 1)*yVec_[idx] +
    (std::pow(t,3) - 2*std::pow(t,2) + t)*(xVec_[idx+1] - xVec_[idx])*mVec_[idx] +
    (-2*std::pow(t,3) + 3*std::pow(t,2))*yVec_[idx+1] +
    (std::pow(t,3) - std::pow(t,2))*(xVec_[idx+1] - xVec_[idx])*mVec_[idx+1];
}


MonotoneCubicInterpolation::MonotoneCubicInterpolation(const double* xVec, const double* yVec, const size_t size):
  CubicHermiteSpline()
{
  std::vector<double> delta(size, 0);
  std::vector<double> mVec(size, 0);
  for(int i=0; i<size-1; ++i)
    delta[i] = (yVec[i+1] - yVec[i])/(xVec[i+1] - xVec[i]);
  for(int i=1; i<size-1; ++i)
    mVec[i] = (delta[i-1] + delta[i])/2.;
  mVec[0] = delta[0]; mVec[size-1] = delta[size-2];
  for(int i=0; i<size-1; ++i)
  {
    if(std::abs(delta[i]) < 1e-30)
      mVec[i] = mVec[i+1] = 0.;
  }
  this->allocate_mem_(xVec, yVec, &mVec[0], size);
}
} // end namespace etc
