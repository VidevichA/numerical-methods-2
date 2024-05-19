import numpy as np
import data


# x[i] * (u[i+1] - 2u[i] + u[i-1]) / h^2 + (u[i+1] - u[i-1]) / (2h) - 3 * u[i] ≈ 2 + 6 * x[i] - x[i]^2 = f[i]
# x[i]/h^2 + 1/2h
# -2 * x[i] / h**2 -3
# x[i]/h^2 - 1/2h


def finite_difference_method():
    h = (data.x_upper - data.x_lower) / (data.N - 1)
    x = np.linspace(data.x_lower, data.x_upper, data.N)

    A = np.zeros((data.N, data.N))
    A[0][0] = -1/h
    A[0][1] = 1/h
    A[-1][-1] = (h + 1) / h
    A[-1][-2] = -1 / h

    b = np.zeros(data.N)
    b[0] = 0.75
    b[-1] = 7

    for i in range(1, data.N - 1):
        b[i] = 2 + 6 * x[i] - x[i]**2
        A[i][i+1] = x[i] / h**2 + 1 / (2 * h)
        A[i][i] = -2 * x[i] / h**2 - 3
        A[i][i-1] = x[i]/h**2 - 1/(2*h)

    y = np.linalg.solve(A, b)
    return x, y
