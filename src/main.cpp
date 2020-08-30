#include <iostream>
#include <iomanip>
#include <fstream>
#include "hirschfye_qmc.hpp"
#include "blas_lapack.hpp"
#include "etc.hpp"
#include "argparse.hpp"
#include "green_tools.hpp"

int main(int argc, char * argv[])
{
  std::vector<pair_t> options, defaults;
  // env; explanation of env
  options.push_back(pair_t("beta", "inverse temperature"));
  options.push_back(pair_t("U", "coulomb interaction"));
  options.push_back(pair_t("mu", "chemical potential"));
  options.push_back(pair_t("D", "half bandwidth"));
  options.push_back(pair_t("V", "hybridization strength"));
  options.push_back(pair_t("nmeas", "# of measurements"));
  options.push_back(pair_t("ntau", "# of slice of imaginary time"));
  options.push_back(pair_t("niwn", "# of Matsubara frequency"));
  options.push_back(pair_t("nwarm", "# of Monte Carlo steps for warming up"));
  options.push_back(pair_t("nms", "# of MCMC steps for sampling Green's function"));
  options.push_back(pair_t("nperiod", "period of clean update"));
  options.push_back(pair_t("seed", "random number seed"));
  options.push_back(pair_t("path", "directory to save a file"));
  // env; default value
  defaults.push_back(pair_t("D", "1.0"));
  defaults.push_back(pair_t("V", "1.0"));
  defaults.push_back(pair_t("ntau", "121"));
  defaults.push_back(pair_t("niwn", "60"));
  defaults.push_back(pair_t("nwarm", "100"));
  defaults.push_back(pair_t("nms", "3"));
  defaults.push_back(pair_t("nperiod", "3"));
  defaults.push_back(pair_t("seed", "1"));
  defaults.push_back(pair_t("path", "."));

  argsparse parser(argc, argv, options, defaults);
  const double beta = parser.find<double>("beta"),
    U = parser.find<double>("U"),
    mu = parser.find<double>("mu"),
    V = parser.find<double>("V"),
    D = parser.find<double>("D");
  const int ntau = parser.find<int>("ntau"),
    niwn = parser.find<int>("niwn"),
    nwarmup = parser.find<int>("nwarm"),
    nms = parser.find<int>("nms"),
    nmeas = parser.find<int>("nmeas"),
    nperiodSweeps = parser.find<int>("nperiod");
  const unsigned long seed = parser.find<unsigned long>("seed");

  // print info of the registered args
  parser.print(std::cout);

  ImTimeGreen gtau(ntau);
  ImFreqGreen giwn = SemiCircularGreen(niwn, beta, mu, V, D);
  GreenInverseFourier(giwn, beta, {1.0, 0.0, 0.0}, gtau);

  // # of mesh (time): ntau -> ntau-1
  ImTimeGreen gtau0(ntau-1), gtau0_up(ntau-1, 0.0), gtau0_dw(ntau-1, 0.0),
    gtau2_up(ntau-1, 0.0), gtau2_dw(ntau-1, 0.0);
  for (int i=0; i<gtau0.size(); ++i)
    gtau0[i] = gtau[i];

  HirschFyeQMC qmc(gtau0.data(), gtau0.size(), beta, U, seed);
  std::cout << "# start Hirsch-Fye QMC" << std::endl;
  qmc.do_montecarlo_step(nwarmup, nperiodSweeps);

  for (int n=0; n<nmeas; ++n)
  {
    std::cout << "\r# loop : " << std::setw(4) << (n+1) << " ";
    const double acceptanceRate = qmc.do_montecarlo_step(nms, nperiodSweeps);
    std::cout << " --- acceptance rate : " << acceptanceRate << std::flush;
    qmc.accumulate(gtau0_up.data(), gtau0_dw.data(), gtau2_up.data(), gtau2_dw.data());
  }
  std::cout << std::endl;

  for (int i=0; i<gtau0.size(); ++i)
  {
    gtau0_up[i] /= nmeas;
    gtau0_dw[i] /= nmeas;
    gtau2_up[i] /= nmeas;
    gtau2_dw[i] /= nmeas;
  }

  // # of mesh (time): ntau-1 -> ntau
  ImTimeGreen gtau_up(ntau), gtau_dw(ntau), chi_up(ntau), chi_dw(ntau);
  for (int i=0; i<ntau-1; ++i)
  {
    gtau_up[i] = -gtau0_up[i];
    gtau_dw[i] = -gtau0_dw[i];
    const double chi1 = std::sqrt((gtau2_up[i]-gtau0_up[i]*gtau0_up[i])/((nmeas-1)));
    const double chi2 = std::sqrt((gtau2_dw[i]-gtau0_dw[i]*gtau0_dw[i])/((nmeas-1)));
    chi_up[i] = chi1;
    chi_dw[i] = chi2;
  }
  gtau_up[ntau-1] = -1.0-gtau_up[0];
  gtau_dw[ntau-1] = -1.0-gtau_dw[0];
  chi_up[ntau-1] = chi_up[0];
  chi_dw[ntau-1] = chi_dw[0];

  std::string Ustr = etc::to_string(U),
    betastr = etc::to_string(beta),
    mustr = etc::to_string(mu);
  const std::string filename1 = parser.find<>("path") + "/Gtau-U" + Ustr + "B" + betastr + "M" + mustr + ".dat";

  std::ofstream file(filename1);
  file << "#       tau        Gtau sigma(Gtau)" << std::endl;
  for (int i=0; i<ntau; ++i)
  {
    const double tau = beta/(ntau-1)*i;
    file << std::setw(11) << tau << " "
	 << std::setw(11) << gtau_up[i] << " "
	 << std::setw(11) << chi_up[i] << " "
	 << std::setw(11) << gtau_dw[i] << " "
	 << std::setw(11) << chi_dw[i] << std::endl;
  }
  file.close();

  return 0;
}
