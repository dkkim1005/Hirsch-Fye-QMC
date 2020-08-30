#include "green_tools.hpp"

void GreenInverseFourier(const ImFreqGreen & giwn, const double beta,
  const std::vector<double> & M, ImTimeGreen & gtau)
{
  // check whether mesh size satisfies the Nyquist theorem.
  assert(gtau.size()/2 >= giwn.size());
  std::vector<std::complex<double> > iwn(giwn.size());
  for (int n=0; n<giwn.size(); ++n)
    iwn[n] = std::complex<double>(0, (2*n+1)*M_PI/beta);
  std::vector<std::complex<double> > giwn_temp(gtau.size(), std::complex<double>(0, 0)),
    gtau_temp(gtau.size());
  // giwn_temp[i] is 0.0 for i >= giwn.size().
  for (int n=0; n<giwn.size(); ++n)
    giwn_temp[n] = giwn[n] - (M[0]/iwn[n] + M[1]/std::pow(iwn[n], 2) + M[2]/std::pow(iwn[n], 3));
  etc::fft_1d_complex(giwn_temp.data(), gtau_temp.data(), gtau.size()-1, etc::FORWARD);	
  for (int i=0; i<gtau.size()-1; ++i)
  {
    const double tau = beta*i/(gtau.size()-1.0);
    gtau[i] = 2.0/beta*(std::exp(std::complex<double>(0.0, -M_PI*tau/beta))*gtau_temp[i]).real()
      - 0.5*M[0] + (tau/2.0-beta/4.0)*M[1] + (tau*beta/4.0 - std::pow(tau, 2)/4.0)*M[2];
  }
  std::complex<double> gedge(0, 0); // := G(beta)
  // giwn_temp[i] is 0.0 for i >= giwn.size().
  for (int n=0; n<giwn.size(); ++n)
    gedge -= giwn_temp[n];
  gedge *= 2.0/beta;
  gtau[gtau.size()-1] = gedge.real() - 0.5*M[0] + (beta/4.0)*M[1];
}

void GreenFourier(const ImTimeGreen & gtau, const double beta, ImFreqGreen & giwn, const bool useSpline)
{
  ImTimeGreen gtautmp = gtau;
  // check whether mesh size satisfies the Nyquist theorem.
  if ((gtau.size()/2 < giwn.size()) or useSpline)
  {
    std::vector<double> tau(gtau.size());
    for (int i=0; i<tau.size(); ++i)
      tau[i] = i*beta/(tau.size()-1.0);
    etc::MonotoneCubicInterpolation interp(tau.data(), gtau.data(), tau.size());
    const int size = 10*giwn.size()+1;
    gtautmp.resize(size);
    for (int i=0; i<gtautmp.size(); ++i)
      gtautmp[i] = interp(i*beta/(gtautmp.size()-1.0));
  }

  const double dtau = beta/(gtautmp.size()-1.0);
  std::vector<std::complex<double> > gmod(gtautmp.size(), 0.0);

  for(int i=0; i<gtautmp.size(); ++i)
  {
    const double tau = i*dtau;
    gmod[i] = gtautmp[i]*std::exp(std::complex<double>(0, M_PI/beta*tau))*dtau;
  }

  // composite Simpson's rule
  for (int i=1; i<gmod.size()-1; i+=2)
    gmod[i] *= 4;
  for (int i=2; i<gmod.size()-2; i+=2)
    gmod[i] *= 2;
  for (int i=0; i<gmod.size(); i++)
    gmod[i] *= 1.0/3.0;

  if (gmod.size()%2 == 0)
  {
    for (int i=gmod.size()-4; i<gmod.size(); ++i)
    {
      const double tau = i*dtau;
      gmod[i] = gtautmp[i]*std::exp(std::complex<double>(0, M_PI/beta*tau))*dtau;
    }
    // Simpson's 3/8 rule
    gmod[gmod.size()-4] *= (1.0/3.0+3.0/8.0);
    gmod[gmod.size()-3] *= (9.0/8.0);
    gmod[gmod.size()-2] *= (9.0/8.0);
    gmod[gmod.size()-1] *= (3.0/8.0);
  }

  std::vector<std::complex<double> > giwn_temp(gtautmp.size(), std::complex<double>(0,0));

  etc::fft_1d_complex(gmod.data(), giwn_temp.data(), gtautmp.size()-1, etc::BACKWARD);
	
  for(int i=0;i<giwn_temp.size()-1;++i)
    giwn_temp[i] += gmod[gmod.size()-1];

  for (int n=0; n<giwn.size(); ++n)
    giwn[n] = giwn_temp[n];
}

ImFreqGreen SemiCircularGreen(const int niwn, const double beta,
  const double mu, const double V, const double D, const int ne)
{
  // hybridization term with semicircular DOS: gamma_n := \rho(e_k)*V^2/(iw_n - e_k)
  struct GammaFunc
  {
    typedef std::complex<double> returnType;
    GammaFunc(const double V, const double D): V_(V), D_(D) {}
    void set_ImFreq(const std::complex<double> iwn) { iwn_ = iwn; }
    // DOS(e) = 2*\sqrt(1-(e/D)^2)/(pi*D)
    returnType operator()(const double & e) const
    {
      return std::pow(V_, 2)*2.0*std::sqrt(1.0-std::pow(e/D_, 2))/(M_PI*D_)/(iwn_ - e);
    }
  private:
    std::complex<double> iwn_;
    const double V_, D_;
  };

  GammaFunc gm(V, D);
  // ek: energy in the conduction band
  std::vector<double> ek(ne);
  for (int i=0; i<ek.size(); ++i)
    ek[i] = -D + i*2.0*D/(ek.size()-1.0);
  ImFreqGreen giwn(niwn);
  for (int n=0; n<niwn; ++n)
  {
    const std::complex<double> iwn = std::complex<double>(0, (2*n+1)*M_PI/beta);
    gm.set_ImFreq(iwn);
    // giwn_n = (iw_n + mu - gamma_n)^-1
    giwn[n] = 1.0/(iwn + mu - etc::simpson_rule(ek, gm));
  }
  return giwn;
};
