#include "hirschfye_qmc.hpp"

HirschFyeQMC::HirschFyeQMC(const double * gtau0, const int ntau, const double beta, const double U, const unsigned long seed):
  nsweeps_(0),
  nchain_(0),
  N_(ntau),
  lambda_(std::log(std::exp(beta/ntau*U/2.0) + std::sqrt(std::pow(std::exp(beta/ntau*U/2.0), 2)-1))),
  gtau0mat_(ntau*ntau, 0.0),
  gtau_up_(ntau*ntau, 0.0),
  gtau_dw_(ntau*ntau, 0.0),
  Gtau_up_(ntau, 0.0),
  Gtau_dw_(ntau, 0.0),
  vi_up_(ntau, 0.0),
  uj_up_(ntau, 0.0),
  vi_dw_(ntau, 0.0),
  uj_dw_(ntau, 0.0),
  e_ls_(ntau, 0.0),
  A_(ntau*ntau, 0.0),
  spins_(ntau, 1)
{
  randdev_.seed(seed);
  // green matrix from vector
  // change a convention of Green's function : G(tau) = -T_t<c(tau)C^+>_\beta --> +T_t<c(tau)C^+>_\beta
  for (int i=0; i<N_; i++)
  {
    for (int j=i; j<N_; j++)
      gtau0mat_[j*N_+i] = gtau0[j-i]*(-1.0);
    //consider antisymmetry (Fermions) : G(beta+tau) = -G(tau)
    for (int j=i+1; j<N_; j++)
      gtau0mat_[i*N_+j] = -gtau0[N_-(j-i)]*(-1.0); 
  }

  for (int i=0; i<spins_.size(); i++)
    spins_[i] = ((unidist_(randdev_)<0.5) ? 1 : -1);

  // memoization of exp calculations
  exp1lamb_[0] = std::exp(-lambda_),
  exp1lamb_[2] = std::exp(lambda_),
  exp2lamb_[0] = std::exp(-2.0*lambda_),
  exp2lamb_[2] = std::exp(2.0*lambda_);

  this->clean_update_(gtau_up_, UP);
  this->clean_update_(gtau_dw_, DOWN);
}

double HirschFyeQMC::do_montecarlo_step(const int nmcsteps, const int nperiodSweeps)
{
  for (int n=0; n<nmcsteps; ++n)
  {
    // do local update
    for (int i=0; i<N_; i++)
    {
      if (this->single_flip_update_())
      {
        nsweeps_ += 1;
        if (nsweeps_ >= nperiodSweeps)
        {
          this->clean_update_(gtau_up_, UP);
          this->clean_update_(gtau_dw_, DOWN);
          nsweeps_ = 0;
        }
      }
    }
  }
  const double acceptanceRateMean = (acceptanceRate_/nchain_);
  acceptanceRate_ = 0.0;
  nchain_ = 0;
  return acceptanceRateMean;
}

void HirschFyeQMC::accumulate(double * G_up, double * G_dw, double * G2_up, double * G2_dw)
{
  this->green_vector_from_matrix_(Gtau_up_, gtau_up_);
  this->green_vector_from_matrix_(Gtau_dw_, gtau_dw_);
  for (int i=0; i<N_; ++i)
  {
    G_up[i] += Gtau_up_[i];
    G_dw[i] += Gtau_dw_[i];
    G2_up[i] += Gtau_up_[i]*Gtau_up_[i];
    G2_dw[i] += Gtau_dw_[i]*Gtau_dw_[i];
  }
}

void HirschFyeQMC::clean_update_(std::vector<double> & gtaumat, SPIN_STATE SPIN)
{
  for(int i=0; i<N_; ++i)
    e_ls_[i] = e1lamb_(SPIN*spins_[i]);
  // calculate A = 1 + (1-g)*(exp(v')-1)
  std::fill(A_.begin(), A_.end(), 0.0);
  for (int i=0; i<N_; i++)
  {
    for (int j=0; j<N_; j++)
      A_[i*N_+j] = -gtau0mat_[i*N_+j]*(e_ls_[j]-1);
    A_[i*N_+i] += e_ls_[i];
  }
  // solve A*gtaumat=gtau0mat
  // gesv solves A*X=B, input for B is gtau0mat, output (=solution X) is gtau
  gtaumat = gtau0mat_;
  etc::transpose(gtaumat.data(), N_);
  etc::transpose(A_.data(), N_);
  lapack::gesv(N_, N_, A_.data(), gtaumat.data());
  etc::transpose(gtaumat.data(), N_);
}

bool HirschFyeQMC::single_flip_update_()
{
  nchain_ += 1;
  // choose random site
  int site = (unsigned int)(N_*unidist_(randdev_));
  // calculate ratio of determinants
  const double r_up = 1 + (1-gtau_up_[site*N_+site])*(e2lamb_(-spins_[site])-1),
    r_down = 1 + (1-gtau_dw_[site*N_+site])*(e2lamb_(spins_[site])-1),
    det_rat = r_up*r_down;

  // heat-bath method: acceptance rate := det_rat/(1+det_rat)
  const double acceptanceRate = std::abs(det_rat/(1+det_rat));
  acceptanceRate_ += acceptanceRate;
  bool isSpinFliped = unidist_(randdev_) < acceptanceRate;

  // if update successful ...
  if (isSpinFliped)
  {
    // update Green's function
    const double tmp1 = e2lamb_(-spins_[site])-1,
      tmp2 = e2lamb_(spins_[site])-1,
      alpha_up = tmp1/(1+(1-gtau_up_[site*N_+site])*tmp1),
      alpha_dw = tmp2/(1+(1-gtau_dw_[site*N_+site])*tmp2);
    for(int i=0; i<N_; ++i)
    {
      vi_up_[i] = gtau_up_[i*N_+site]-(i==(int)site);
      vi_dw_[i] = gtau_dw_[i*N_+site]-(i==(int)site);
      uj_up_[i] = gtau_up_[site*N_+i];
      uj_dw_[i] = gtau_dw_[site*N_+i];
    }

    blas::ger(N_, N_, alpha_up, vi_up_.data(), uj_up_.data(), gtau_up_.data());
    blas::ger(N_, N_, alpha_dw, vi_dw_.data(), uj_dw_.data(), gtau_dw_.data());

    // update spin
    spins_[site] = -spins_[site];
  }

  return isSpinFliped;
}

void HirschFyeQMC::green_vector_from_matrix_(std::vector<double> & gtauvec, const std::vector<double> & gtaumat)
{
  std::fill(gtauvec.begin(), gtauvec.end(), 0.0);
  for (int i=0; i<N_; i++)
  {
    for (int j=i; j<N_; j++)
    {
      gtauvec[j-i] += gtaumat[j*N_+i];
      if (j>i) //consider antisymmetry (Fermions) : G(beta+tau) = -G(tau)
        gtauvec[N_-(j-i)] -= gtaumat[i*N_+j];
    }
  }
  for(int i=0; i<N_; ++i)
    gtauvec[i] /= N_;
}
