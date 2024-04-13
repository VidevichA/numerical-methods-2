import numpy as np


def runge_kutta_third_order(f, h, y0, x_values):
    y_values = np.zeros((len(x_values), len(y0)))
    y_values[0] = y0
    for i in range(1, len(x_values)):
        xi = x_values[i]
        yi_prev = y_values[i - 1]
        fi0 = f(xi, yi_prev)
        fi1 = f(xi + h / 3, yi_prev + h / 3 * fi0)
        fi2 = f(xi + 2 * h / 3, yi_prev + 2 * h / 3 * fi1)
        yi_new = yi_prev + h / 4 * (fi0 + 3 * fi2)
        y_values[i] = yi_new
    return x_values, y_values
