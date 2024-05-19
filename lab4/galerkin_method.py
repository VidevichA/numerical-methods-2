import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import data


def fi(k, x):
    if k == 0:
        return 0.75 * x + 5.5
    else:
        return -(data.x_upper-data.x_lower)**(k+1)-(k+1)*(data.x_upper-data.x_lower)**k+(x-data.x_lower)**(k+1)


def d_fi(k, x):
    h = 0.001
    return (fi(k, x+h)-fi(k, x-h))/(2*h)


def d2_fi(k, x):
    h = 0.001
    return (d_fi(k, x+h)-d_fi(k, x-h))/(2*h)


def r(x):
    return 1/x


def q(x):
    return -3/x


def integrate(f):
    return spi.quad(f, data.x_lower, data.x_upper)[0]


def get_func_to_integrate_fi(i):
    def func(x):
        return (d2_fi(0, x)+r(x)*d_fi(0, x)+q(x)*fi(0, x) - 2/x - 6 + x)*fi(i, x)
    return func


def f(i):
    return -1 * integrate(get_func_to_integrate_fi(i))


def get_func_to_integrate_Aij(i, j):
    def func(x):
        return (d2_fi(j, x)+r(x)*d_fi(j, x)+q(x)*fi(j, x))*fi(i, x)
    return func


def galerkin_method():
    n = 10
    G = []
    for i in range(1, n):
        tmp = []
        for j in range(1, n):
            tmp.append(integrate(get_func_to_integrate_Aij(i, j)))
        G.append(tmp[:])

    K = []

    for i in range(1, n):
        K.append(f(i))

    G = np.asarray(G, dtype=np.float64)
    K = np.asarray(K, dtype=np.float64)

    coeff = np.linalg.solve(G, K)

    def u(x):
        sum = fi(0, x)
        for i in range(n-1):
            sum += coeff[i]*fi(i+1, x)
        return sum

    x = np.linspace(data.x_lower, data.x_upper, data.N)
    y = map(lambda x: u(x), x)
    return x, list(y)
