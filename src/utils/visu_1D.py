import matplotlib.pyplot as plt


def plot_solution(grid, approximate_sol, exact_sol, title="Heat eq"):
    """

    :param grid:
    :param approximate_sol:
    :param exact_sol:
    :param title:
    :return:
    """
    plt.plot(grid, exact_sol, label="exact solution")
    plt.plot(grid, approximate_sol, label="approximate solution")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.title(title)
    plt.legend()
    plt.show()
