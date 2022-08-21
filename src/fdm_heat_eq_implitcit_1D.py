import math
import os
import numpy

from src.utils import visu1D
from ctypes import POINTER, c_int, c_float, cdll

"""
MacOS: Before running this script, please build fortran library [1]:
    gfortran -dynamiclib src/utils/tridiag.f90 -o src/tridiag.dylib
There will be a file called 'tridiag.dylib' in the project src folder.
Reference:
[1] http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tridiag_f90.txt
[2] https://stackoverflow.com/questions/19263879/speeding-up-element-wise-array-multiplication-in-python/19458585#19458585
"""
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
FORTRAN_LIB_PATH = os.path.join(ROOT_DIR, 'tridiag.dylib')
fortran = cdll.LoadLibrary(FORTRAN_LIB_PATH)
fortran.tridiag.argtypes = [POINTER(c_float),
                            POINTER(c_float),
                            POINTER(c_float),
                            POINTER(c_float),
                            POINTER(c_float),
                            POINTER(c_int),
                            POINTER(c_int)]
fortran.tridiag.restype = None

if __name__ == "__main__":
    # FDM, implicit scheme for the heat equation:
    # INPUT
    # Advection velocity:
    mu = 1 / 16
    # Grid size; Use periodic boundary conditions
    X, M, Tend = 1, 10, 0.5
    h = 1 / (M + 1)
    Dt = 0.02
    x = numpy.linspace(0, X, M + 2)  # x(1) =0, x(2) = h, , ..., x(M) = X -h; x(M+1) = X;
    wn = numpy.sin(2 * math.pi * x)
    wm = numpy.sin(2 * math.pi * x)
    time = 0
    mu = mu * Dt / h ** 2

    A = numpy.zeros(M + 2)
    A[0] = 1.0 + 2.0 * mu
    A[1:] = -1.0 * mu

    C = numpy.zeros(M + 2)
    C[0:M + 1] = -1.0 * mu
    C[-1] = 1.0 + 2.0 * mu

    B = numpy.zeros(M + 2)
    B[:] = 1.0 + 2.0 * mu

    while time < Tend:
        a, b, c = A, B, C
        code = -1
        fortran.tridiag(a.ctypes.data_as(POINTER(c_float)),
                        b.ctypes.data_as(POINTER(c_float)),
                        b.ctypes.data_as(POINTER(c_float)),
                        wm.ctypes.data_as(POINTER(c_float)),
                        wn.ctypes.data_as(POINTER(c_float)),
                        c_int(M + 2), c_int(code))
        wn = wm

        time += Dt

    w = numpy.exp(-time / 4 * math.pi ** 2) * numpy.sin(2 * math.pi * x)
    visu1D.plot_solution(x, wn, w)
