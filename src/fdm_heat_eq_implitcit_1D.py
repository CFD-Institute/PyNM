import math
import numpy

from src.utils import visu1D


def thomas(
        a: numpy.ndarray,
        b: numpy.ndarray,
        c: numpy.ndarray,
        r: numpy.ndarray,
        u: numpy.ndarray,
        n: int):
    """

    :param a:
    :param b:
    :param c:
    :param r:
    :param u:
    :param n:
    :return:
    """
    eps = 10**(-9)

    if b[0] <= eps:
        return

    bet = b[0]
    u[0] = r[0] / bet

    gam = numpy.zeros(n)

    for j in range(1, n):
        gam[j] = c[j - 1] / bet
        bet = b[j] - a[j] * gam[j]
        if b[0] <= eps:
            return
        u[j] = (r[j] - a[j] * u[j - 1]) / bet

    for j in range(n - 2, -1, -1):
        u[j] = u[j] - gam[j + 1] * u[j + 1]


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

    while time < Tend:
        a, b, c = A, B, C
        thomas(a, b, c, wm, wn, M + 2)
        wn = wm

        time += Dt

    w = numpy.exp(-time / 4 * math.pi ** 2) * numpy.sin(2 * math.pi * x)
    visu1D.plot_solution(x, wm, w)
