import math
import numpy
import time as timer

from src.libs.tridiag import TriDiagonalMatrix
from src.utils import visu1D


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

    start = timer.time()
    while time < Tend:
        a, b, c = A, B, C
        TriDiagonalMatrix.thomas(a, b, c, wm, wn, M + 2)
        wn = wm

        time += Dt

    end = timer.time()
    print("Total elapsed time: ", end - start)

    w = numpy.exp(-time / 4 * math.pi ** 2) * numpy.sin(2 * math.pi * x)
    visu1D.plot_solution(x, wm, w)
