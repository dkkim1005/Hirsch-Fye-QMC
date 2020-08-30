import numpy as np
from . import _cpp_module


class hfqmc:
    def __init__(self, gtau0, beta, U, seed = 0):
        """
          gtau0 : imaginary time Green's function (np.ndarray, float64)
           beta : inverse temperature (float64)
              U : onsite interaction (float64)
           seed : random number seed (int)
        """
        if type(gtau0) is np.ndarray:
            if gtau0.dtype is not np.dtype('float64'):
                raise TypeError("gtau0.dtype is not dtype(float64)")
        else:
            raise TypeError("type(gtau0) should be 'np.ndarray' type")
        gtau1 = np.array([gtau0[i] for i in range(len(gtau0)-1)])
        self._impl = _cpp_module.hfqmc(gtau1, beta, U, seed)
        self._ntau = len(gtau1)
        self._spins = np.zeros([self._ntau], dtype = 'int')


    def do_montecarlo_step(self, nmcsteps):
        """
          nmcsteps : # of Monte Carlo steps (int)
        """
        if type(nmcsteps) is not int:
            raise TypeError("type(nmcsteps) should be 'int' type")
        self._impl.do_montecarlo_step(nmcsteps, 3)


    def meas(self, nmcsteps, nmeas):
        """
          nmcsteps : # of Monte Carlo steps (int)
          nmeas : # of measurements (int)
        """
        _gtau1_up = np.zeros([self._ntau], dtype = 'float64')
        _gtau1_dw = np.zeros([self._ntau], dtype = 'float64')
        _gtau2_up = np.zeros([self._ntau], dtype = 'float64')
        _gtau2_dw = np.zeros([self._ntau], dtype = 'float64')
        for n in range(nmeas):
            self._impl.do_montecarlo_step(nmcsteps, 3)
            self._impl.accumulate(_gtau1_up, _gtau1_dw, _gtau2_up, _gtau2_dw)

        _gtau1_up /= -1.0*nmeas
        _gtau1_dw /= -1.0*nmeas
        _gtau2_up /= nmeas
        _gtau2_dw /= nmeas

        _err_up = np.sqrt((_gtau2_up - _gtau1_up**2)/(nmeas-1))
        _err_dw = np.sqrt((_gtau2_dw - _gtau1_dw**2)/(nmeas-1))

        gtau_up = np.zeros([self._ntau+1])
        gtau_up[:self._ntau] = _gtau1_up
        gtau_dw = np.zeros([self._ntau+1])
        gtau_dw[:self._ntau] = _gtau1_dw
        err_up = np.zeros([self._ntau+1])
        err_up[:self._ntau] = _err_up
        err_dw = np.zeros([self._ntau+1])
        err_dw[:self._ntau] = _err_dw
        gtau_up[-1] = -1 - gtau_up[0]
        gtau_dw[-1] = -1 - gtau_dw[0]
        err_up[-1] = err_up[0]
        err_dw[-1] = err_dw[0]

        return gtau_up, err_up, gtau_dw, err_dw


    def get_auxiliary_fields(self):
        self._impl.get_auxiliary_fields(self._spins)
        return self._spins
