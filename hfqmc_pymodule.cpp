#include "hirschfye_qmc.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

typedef py::array_t<double> PyImTimeGreen;

class PythonHirschFyeQMC
{
public:
  PythonHirschFyeQMC(const PyImTimeGreen pygtau0, const double beta, const double U, const long seed):
    hfqmc_(nullptr)
  {
    const auto temp = pygtau0.unchecked<1>();
    const double * ptr = &temp(0);
    hfqmc_ = new HirschFyeQMC(ptr, temp.shape(0), beta, U, seed);
  }

  double do_montecarlo_step(const int nmcsteps, const int nperiodSweeps)
  {
    return hfqmc_->do_montecarlo_step(nmcsteps, nperiodSweeps);
  }

  void accumulate(PyImTimeGreen pyGup, PyImTimeGreen pyGdw, PyImTimeGreen pyG2up, PyImTimeGreen pyG2dw)
  {
    auto Guptemp = pyGup.mutable_unchecked<1>();
    auto Gdwtemp = pyGdw.mutable_unchecked<1>();
    auto G2uptemp = pyG2up.mutable_unchecked<1>();
    auto G2dwtemp = pyG2dw.mutable_unchecked<1>();
    hfqmc_->accumulate(&Guptemp(0), &Gdwtemp(0), &G2uptemp(0), &G2dwtemp(0));
  }

  void get_auxiliary_fields(py::array_t<int> pyspins) const
  {
    auto temp = pyspins.mutable_unchecked<1>();
    int * ptr = &temp(0);
    hfqmc_->get_auxiliary_fields(ptr);
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

PYBIND11_MODULE(cpp_module, m)
{
  m.doc() = "Hirsch-Fye QMC";
  py::class_<PythonHirschFyeQMC>(m, "hfqmc")
    .def(py::init<const PyImTimeGreen, const double, const double, const long>())
    .def("do_montecarlo_step", &PythonHirschFyeQMC::do_montecarlo_step)
    .def("accumulate", &PythonHirschFyeQMC::accumulate)
    .def("get_auxiliary_fields", &PythonHirschFyeQMC::get_auxiliary_fields);
}
