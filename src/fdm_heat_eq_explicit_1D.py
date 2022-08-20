import math
import numpy
import matplotlib.pyplot as plt

"""
Original version was in Matlab. See [1] for further details.
Translate to Python by TgDg.
Reference: [1] Boualem Khouider, Finite difference and finite volume methods for transport and conservation laws.
"""


def plot_solution(x, wn):
    """

    :param x:
    :param wn:
    :return:
    """
    xx = numpy.linspace(0, 1, 1000)
    plt.plot(xx, numpy.exp(-time / 4 * math.pi ** 2) * numpy.sin(2 * math.pi * xx), label="exact solution")
    plt.plot(x, wn, label="approximate solution")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.title("Heat eq")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # FDM, forward difference, explicit scheme for the heat equation:
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

    plot_solution(x, wn)

