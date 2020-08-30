import numpy as np
from . import _cpp_module


def SemiCircularGreen(niwn, beta, mu, V = 1.0, D = 1.0):
    """
      niwn : # of Matsubara frequencies (int)
      beta : inverse temperature (float)
        mu : chemical potential (float)
         V : hybridization strength (float, default : 1.0)
         D : half-bandwidth (float, default : 1.0)
    """
    if type(niwn) is not int:
        raise TypeError("type(niwn) should be 'int' type.")
    return _cpp_module.SemiCircularGreen(niwn, beta, mu, V, D)


def GreenFourier(gtau, niwn, beta):
    """
      gtau : imaginary time Green's function (np.ndarray, float64)
      niwn : # of Matsubara frequencies (int)
      beta : inverse temperature (float)
    """
    if type(gtau) is np.ndarray:
        if gtau.dtype is not np.dtype('float64'):
            raise TypeError("gtau.dtype is not dtype(float64)")
    else:
        raise TypeError("type(gtau) should be 'np.ndarray' type")
    if type(ntau) is not int:
        raise TypeError("type(ntau) should be 'int' type")
    giwn = np.zeros([niwn], dtype = 'complex128')
    _cpp_module.GreenFourier(gtau, beta, giwn)
    return giwn


def GreenInverseFourier(giwn, ntau, beta, M = [1.0, 0.0, 0.0]):
    """
      giwn : Matsubara Green's function (np.ndarray, complex128)
      ntau : # of meshes of imaginary time (int)
      beta : inverse temperature (float)
         M : high-frequency correction M[0]/iwn + M[1]/iwn^2 + M[2]/iwn^3
    """
    if type(giwn) is np.ndarray:
        if giwn.dtype is not np.dtype('complex128'):
            raise TypeError("giwn.dtype is not dtype(complex128)")
    else:
        raise TypeError("type(giwn) should be 'np.ndarray' type")
    if type(ntau) is not int:
        raise TypeError("type(ntau) should be 'int' type")
    gtau = np.zeros([ntau], dtype = np.float64)
    m = np.zeros([3], dtype = 'float64')
    try:
        for i in range(3):
            m[i] = M[i]
    except IndexError:
        pass
    _cpp_module.GreenInverseFourier(giwn, beta, m, gtau)
    return gtau
