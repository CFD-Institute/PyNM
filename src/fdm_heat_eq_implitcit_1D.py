import math
import numpy

from src.utils import visu1D


if __name__ == "__main__":
    # FDM, implicit scheme for the heat equation:
    # INPUT
    # Advection velocity:
    mu = 1 / 16
    # Grid size; Use periodic boundary conditions
    X, M, Tend = 1, 10, 1
    h = 1 / (M + 1)
    Dt = 0.02
    x = numpy.linspace(0, X, M + 2)  # x(1) =0, x(2) = h, , ..., x(M) = X -h; x(M+1) = X;
    wn = numpy.sin(2 * math.pi * x)
    time = 0
    mu = mu * Dt / h ** 2

    while time < Tend:
        wn[1:M + 1] = wn[1:M + 1] + mu * (wn[2:M + 2] - 2 * wn[1:M + 1] + wn[0:M])
        time += Dt

    w = numpy.exp(-time / 4 * math.pi ** 2) * numpy.sin(2 * math.pi * x)
    visu1D.plot_solution(x, wn, w)

