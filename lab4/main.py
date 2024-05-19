import matplotlib.pyplot as plt
from finite_difference_method import finite_difference_method
from galerkin_method import galerkin_method
import numpy as np

x1, u1 = finite_difference_method()
x2, u2 = galerkin_method()

print(np.linalg.norm(np.array(u1) - np.array(u2)))

plt.plot(x1, u1, label='Finite Difference Method')
plt.plot(x2, u2, label='Galerkin Method')

plt.legend()
plt.show()
