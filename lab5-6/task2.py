from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

n = 2
m = 1
h = 0.1
x1 = np.linspace(0, m, int(m/h)+1)
x2 = np.linspace(0, n, int(n/h)+1)

u = np.zeros((len(x1), len(x2)))
fi = np.zeros((len(x1), len(x2)))


def f(x1, x2):
    return -(x2**2 + 2) * np.exp(-x1)


def f1(x1):
    return 0


def f2(x2):
    return x2**2


def f3(x1):
    return 4*np.exp(-x1)


def f4(x2):
    return x2**2*np.exp(-1)


for i in range(len(x1)):
    for j in range(len(x2)):
        fi[i][j] = f(x1[i], x2[j])


for i in range(len(x1)):
    u[i][0] = f1(x1[i])


for i in range(len(x2)):
    u[0][i] = f2(x2[i])


for i in range(len(x1)):
    u[i][-1] = f3(x1[i])


for i in range(len(x2)):
    u[-1][i] = f4(x2[i])

eps = 0.001
u_prev_prev = np.copy(u)
u_prev = np.copy(u)
k = 0
while True:
    for i in range(1, len(x1)-1):
        for j in range(1, len(x2)-1):
            u[i][j] = 1/4 * (u[i-1][j] + u[i+1][j] + u[i]
                             [j-1] + u[i][j+1] + h * h * fi[i][j])
    k += 1
    if k != 1:
        v = np.max(np.abs(u-u_prev)) / np.max(np.abs(u_prev-u_prev_prev))
        if (np.max(np.abs(u-u_prev)) <= eps * (1-v)):
            break
    u_prev_prev = np.copy(u_prev)
    u_prev = np.copy(u)


def plot_3d_matrix(u):
    x = np.linspace(0, 1, u.shape[1])
    y = np.linspace(0, 2, u.shape[0])
    X, Y = np.meshgrid(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, u, cmap='viridis')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()


print(u)
plot_3d_matrix(u)


def real_solution(x1, x2):
    return np.exp(-x1) * x2**2


sol = np.zeros((len(x1), len(x2)))
for i in range(len(x1)):
    for j in range(len(x2)):
        sol[i][j] = real_solution(x1[i], x2[j])


print(np.max(np.abs(sol-u)))
