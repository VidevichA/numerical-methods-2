import numpy as np
from scipy.optimize import fsolve


def implicit_trapezoidal_method(f, h, y0, x_values):
    y_values = np.zeros((len(x_values), len(y0)))
    y_values[0] = y0
    for i in range(1, len(x_values)):
        xi = x_values[i]
        yi_prev = y_values[i - 1]

        def equation(y_next):
            return yi_prev + h / 2 * (f(xi, yi_prev) + f(xi + h, y_next)) - y_next
        yi_next = fsolve(equation, yi_prev)
        y_values[i] = yi_next
    return x_values, y_values
