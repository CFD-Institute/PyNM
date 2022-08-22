import math
import os
import numpy
import time as timer

from src.utils import visu1D
from ctypes import POINTER, c_int, c_double, cdll

"""
MacOS: Before running this script, please build fortran library [1]:
    src$ gfortran -shared -fPIC libs/tridiag.f90 -o tridiag.dylib
There will be a file called 'tridiag.dylib' in the project src folder.
Reference:
[1] http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tridiag_f90.txt
[2] https://stackoverflow.com/questions/19263879/speeding-up-element-wise-array-multiplication-in-python/19458585#19458585
"""
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
FORTRAN_LIB_PATH = os.path.join(ROOT_DIR, 'tridiag.dylib')
fort_lib = cdll.LoadLibrary(FORTRAN_LIB_PATH)
tri_diag = fort_lib.__getattr__("tridiag")
tri_diag.argtypes = [POINTER(c_double),
                     POINTER(c_double),
                     POINTER(c_double),
                     POINTER(c_double),
                     POINTER(c_double),
                     POINTER(c_int),
                     POINTER(c_int)]
tri_diag.restype = None

if __name__ == "__main__":
    # FDM, implicit scheme for the heat equation:
    # INPUT
    # Advection velocity:
    mu = 1 / 16
    # Grid size; Use periodic boundary conditions
    X, M, Tend = 1, 5000, 1
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

    code = c_int(-1)
    n = c_int(M + 2)

    start = timer.time()
    while time < Tend:
        a, b, c = A, B, C
        tri_diag(a.ctypes.data_as(POINTER(c_double)),
                 b.ctypes.data_as(POINTER(c_double)),
                 c.ctypes.data_as(POINTER(c_double)),
                 wm.ctypes.data_as(POINTER(c_double)),
                 wn.ctypes.data_as(POINTER(c_double)),
                 n, code)
        wn = wm

        time += Dt

    end = timer.time()
    print("Total elapsed time: ", end - start)

    w = numpy.exp(-time / 4 * math.pi ** 2) * numpy.sin(2 * math.pi * x)
    visu1D.plot_solution(x, wn, w)
