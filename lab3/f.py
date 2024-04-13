import numpy as np


def f(x, y):
    u, v, w = y
    return np.array([u + v + np.exp(x), u + x * np.exp(x), w + np.exp(x)])
