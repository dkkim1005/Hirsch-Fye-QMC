#include "hirschfye_qmc.hpp"
#include "green_tools.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

typedef py::array_t<double> PyImTimeGreen;
typedef py::array_t<std::complex<double> > PyImFreqGreen;

class PythonHirschFyeQMC
{
public:
  PythonHirschFyeQMC(const PyImTimeGreen pygtau0, const double beta, const double U, const long seed):
    hfqmc_(nullptr)
  {
    const auto wrap_gtau0 = pygtau0.unchecked<1>();
    hfqmc_ = new HirschFyeQMC(&wrap_gtau0(0), wrap_gtau0.shape(0), beta, U, seed);
  }

  double do_montecarlo_step(const int nmcsteps, const int nperiodSweeps)
  {
    return hfqmc_->do_montecarlo_step(nmcsteps, nperiodSweeps);
  }

  void accumulate(PyImTimeGreen pyGup, PyImTimeGreen pyGdw, PyImTimeGreen pyG2up, PyImTimeGreen pyG2dw)
  {
    auto wrap_Gup = pyGup.mutable_unchecked<1>();
    auto wrap_Gdw = pyGdw.mutable_unchecked<1>();
    auto wrap_G2up = pyG2up.mutable_unchecked<1>();
    auto wrap_G2dw = pyG2dw.mutable_unchecked<1>();
    hfqmc_->accumulate(&wrap_Gup(0), &wrap_Gdw(0), &wrap_G2up(0), &wrap_G2dw(0));
  }

  void get_auxiliary_fields(py::array_t<int> pyspins) const
  {
    auto wrap_spins = pyspins.mutable_unchecked<1>();
    hfqmc_->get_auxiliary_fields(&wrap_spins(0));
  }

  ~PythonHirschFyeQMC()
  {
    if (hfqmc_ != nullptr)
      delete hfqmc_;
    hfqmc_ = nullptr;
  }

private:
  HirschFyeQMC * hfqmc_;
};


void PyGreenInverseFourier(const PyImFreqGreen pygiwn, const double beta,
  const py::array_t<double> pyM, PyImTimeGreen pygtau)
{
  auto wrap_giwn = pygiwn.unchecked<1>();
  auto wrap_M = pyM.unchecked<1>();
  auto wrap_gtau = pygtau.mutable_unchecked<1>();
  ImFreqGreen giwn(&wrap_giwn(0), &wrap_giwn(0) + wrap_giwn.shape(0));
  std::vector<double> M(&wrap_M(0), &wrap_M(0) + wrap_M.shape(0));
  ImTimeGreen gtau(&wrap_gtau(0), &wrap_gtau(0) + wrap_gtau.shape(0));
  GreenInverseFourier(giwn, beta, M, gtau);
  for (int i=0; i<gtau.size(); ++i)
    wrap_gtau(i) = gtau[i];
}

void PyGreenFourier(const PyImTimeGreen pygtau, const double beta, PyImFreqGreen pygiwn)
{
  auto wrap_gtau = pygtau.unchecked<1>();
  auto wrap_giwn = pygiwn.mutable_unchecked<1>();
  ImTimeGreen gtau(&wrap_gtau(0), &wrap_gtau(0) + wrap_gtau.shape(0));
  ImFreqGreen giwn(&wrap_giwn(0), &wrap_giwn(0) + wrap_giwn.shape(0));
  GreenFourier(gtau, beta, giwn, true);
  for (int n=0; n<giwn.size(); ++n)
    wrap_giwn[n] = giwn[n];
}

PyImFreqGreen PySemiCircularGreen(const int niwn, const double beta,
  const double mu, const double V, const double D)
{
  const ImFreqGreen giwn = SemiCircularGreen(niwn, beta, mu, V, D);
  return py::array_t<std::complex<double> >(giwn.size(), giwn.data());
}

PYBIND11_MODULE(_cpp_module, m)
{
  m.doc() = "c++ implementation of Hirsch-Fye QMC";
  py::class_<PythonHirschFyeQMC>(m, "hfqmc")
    .def(py::init<const PyImTimeGreen, const double, const double, const long>())
    .def("do_montecarlo_step", &PythonHirschFyeQMC::do_montecarlo_step)
    .def("accumulate", &PythonHirschFyeQMC::accumulate)
    .def("get_auxiliary_fields", &PythonHirschFyeQMC::get_auxiliary_fields);
  m.def("GreenInverseFourier", &PyGreenInverseFourier, "inverse fourier transform of a Green's function");
  m.def("GreenFourier", &PyGreenFourier, "fourier transform of a Green's function");
  m.def("SemiCircularGreen", &PySemiCircularGreen, "bare Green's function with semi circular DOS");
}
