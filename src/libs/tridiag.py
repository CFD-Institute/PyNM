import numpy


class TriDiagonalMatrix:

    @classmethod
    def thomas(
            cls,
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