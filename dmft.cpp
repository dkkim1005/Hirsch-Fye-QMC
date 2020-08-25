#include <iostream>
#include <iomanip>
#include <fstream>
#include "hirschfye_qmc.hpp"
#include "blas_lapack.hpp"
#include "etc.hpp"
#include "argparse.hpp"
#include "green_tools.hpp"

double distance(const ImTimeGreen & G1, const ImTimeGreen & G2)
{
  double tmp = 0;
  for (int i=0; i<G1.size(); ++i)
    tmp += std::pow(G1[i]-G2[i], 2);
  return std::sqrt(tmp/G1.size());
}

int main(int argc, char * argv[])
{
  std::vector<pair_t> options, defaults;
  // env; explanation of env
  options.push_back(pair_t("beta", "inverse temperature"));
  options.push_back(pair_t("U", "coulomb interaction"));
  options.push_back(pair_t("mu", "chemical potential"));
  options.push_back(pair_t("D", "half bandwidth"));
  options.push_back(pair_t("nloop", "# of dmft loops"));
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
    D = parser.find<double>("D");
  const int ntau = parser.find<int>("ntau"),
    niwn = parser.find<int>("niwn"),
    nwarmup = parser.find<int>("nwarm"),
    nms = parser.find<int>("nms"),
    nloop = parser.find<int>("nloop"),
    nmeas = parser.find<int>("nmeas"),
    nperiodSweeps = parser.find<int>("nperiod");
  const unsigned long seed = parser.find<unsigned long>("seed");

  // print info of the registered args
  parser.print(std::cout);

  ImTimeGreen bare_gtau(ntau), full_gtau(ntau), full_gtau1(ntau);
  ImFreqGreen bare_giwn(niwn), full_giwn(niwn);

  full_giwn = SemiCircularGreen(niwn, beta, 0.0, 1.0, D);
  GreenInverseFourier(full_giwn, beta, {1.0, 0.0, 0.0}, full_gtau1);

  // # of mesh (time): ntau -> ntau-1
  ImTimeGreen bare_gtau0(ntau-1), full_gtau0_up(ntau-1, 0.0), full_gtau0_dw(ntau-1, 0.0),
    full_gtau02_up(ntau-1, 0.0), full_gtau02_dw(ntau-1, 0.0);

  for (int nl=0; nl<nloop; nl++)
  {
    std::cout << "# self-consistent loop: " << nl << std::endl;
    // g(iwn) = 1/(iwn + mu - \frac{D^2}{4}G(iwn))
    for (int n=0; n<niwn; ++n)
      bare_giwn[n] = 1.0/(std::complex<double>(0.0, (2*n+1)*M_PI/beta) + mu - std::pow(D/2.0, 2)*full_giwn[n]);
    GreenInverseFourier(bare_giwn, beta, {1.0, 0.0, 0.0}, bare_gtau);

    for (int i=0; i<bare_gtau0.size(); ++i)
      bare_gtau0[i] = bare_gtau[i];

    HirschFyeQMC qmc(bare_gtau0, beta, U, seed);
    qmc.do_montecarlo_step(nwarmup, nperiodSweeps);

    for (int n=0; n<nmeas; ++n)
    {
      std::cout << "\r# qmc loop : " << std::setw(4) << (n+1) << " " << nmeas;
      const double acceptanceRate = qmc.do_montecarlo_step(nms, nperiodSweeps);
      std::cout << " --- acceptance rate : " << acceptanceRate << std::flush;
      qmc.accumulate(full_gtau0_up, full_gtau0_dw, full_gtau02_up, full_gtau02_dw);
    }
    std::cout << std::endl;

    for (int i=0; i<bare_gtau0.size(); ++i)
    {
      full_gtau0_up[i] /= nmeas;
      full_gtau0_dw[i] /= nmeas;
      full_gtau02_up[i] /= nmeas;
      full_gtau02_dw[i] /= nmeas;
    }

    for (int i=0; i<ntau-1; ++i)
      full_gtau[i] = -(full_gtau0_up[i] + full_gtau0_dw[i])/2.0;
    full_gtau[ntau-1] = -1.0-full_gtau[0];

    const double dist = distance(full_gtau, full_gtau1);
    std::cout << "# dist : " << dist << std::endl;
    if (distance(full_gtau, full_gtau1) < 1e-3)
    {
      std::cout << "# converge!" << std::endl;
      break;
    }
    full_gtau1 = full_gtau;

    GreenFourier(full_gtau, beta, full_giwn);
  }

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
	 << std::setw(11) << full_gtau[i] << std::endl;
  }
  file.close();

  return 0;
}
