import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

b = 1.0
bt = 0.5
N = 20
M = 400
h = b / N
r = bt / M


x = np.linspace(0, b, N+1)
t = np.linspace(0, bt, M+1)

u = np.zeros((N+1, M+1))


def f(x, t):
    return -np.sin(6*x*t)


u[:, 0] = np.ones(N+1)
u[0, :] = np.ones(M+1)


def g1(t):
    return 1


def g2(t):
    return np.cos(t)


for n in range(0, M):
    for i in range(1, N):
        u[i, n+1] = u[i, n] + r / h**2 * \
            (u[i+1, n] - 2*u[i, n] + u[i-1, n]) + r * f(x[i], t[n])
    u[N, n+1] = u[N-1, n+1] + h * g2(t[n+1])


X, T = np.meshgrid(x, t)
U = u.T

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, T, U, cmap='viridis')

ax.set_xlabel('X (пространственная координата)')
ax.set_ylabel('T (время)')
ax.set_zlabel('Температура u(x,t)')
ax.set_title('Распределение температуры в стержне')
plt.show()
