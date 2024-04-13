import numpy as np


def implicit_trapezoidal_method(f, h, y0, x_values):
    y_values = np.zeros((len(x_values), len(y0)))
    y_values[0] = y0
    for i in range(1, len(x_values)):
        xi = x_values[i]
        yi = y_values[i - 1]
        yi_temp = yi + h * 0.5 * \
            (f(xi, yi) + f(xi + h, yi + h * f(xi, yi)))
        y_values[i] = yi_temp
    return x_values, y_values
