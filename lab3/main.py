from implicit_trapezoidal_method import implicit_trapezoidal_method
from runge_kutta_third_order import runge_kutta_third_order
from input_data import *
from f import f
from display_data import display_data

x_values, y_values = implicit_trapezoidal_method(f, h, u0, x_values)
display_data(x_values, y_values,
             'Решение системы дифференциальных уравнений неявным методом трапеций')

y_imt_right = y_values[-1]

x_values, y_values = runge_kutta_third_order(f, h, u0, x_values)
display_data(x_values, y_values,
             'Решение системы дифференциальных уравнений методом Рунге-Кутты 3-го порядка')

y_rk_right = y_values[-1]

print("Модуль разности решений в крайней правой точке интервала: ",
      abs(y_imt_right-y_rk_right))
