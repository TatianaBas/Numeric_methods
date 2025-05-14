"""
Метод стрельбы:
 y(2) = 1.5    y(2) ≈ mu = 1.5000018477666897
Ошибка Рунге–Ромберга: 0.14286168952199954

Метод конечных разностей:
Ошибка Рунге–Ромберга: 2.155419955252036e-05
"""

import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

# —————————————————————— График ——————————————————————
def Plot2(Xi, Yi, X, Y, title=''):
    plt.figure(figsize=(8, 6))
    
    # Численные точки
    plt.scatter(Xi, Yi, c='red', label='Численное решение', zorder=3)

    # Нумерация шагов
    for i, (x, y) in enumerate(zip(Xi, Yi)):
        plt.text(x, y + 0.03, f'{i}', fontsize=8, color='darkred', ha='center')

    # Точное решение
    plt.plot(X, Y, c='blue', label='Точное решение', linewidth=2)

    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def PlotShootingIterations(all_paths, X_exact, Y_exact):
    plt.figure(figsize=(8, 6))
    
    for i, (Xk, Yk) in enumerate(all_paths):
        plt.plot(Xk, Yk, label=f'Итерация {i}', linewidth=1.5, linestyle='--')

    plt.plot(X_exact, Y_exact, c='blue', label='Точное решение', linewidth=2)
    plt.title("Метод стрельбы – итерации")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# —————————————————————— Метод Рунге–Ромберга ——————————————————————
def RungeRombergMet(Met, h, r, p):
    y1 = Met.Get(h)
    y2 = Met.Get(h * r)
    y2_reduced = y2[::r]

    min_len = min(len(y1), len(y2_reduced))
    max_error = 0
    for i in range(min_len):
        err = abs((y1[i] - y2_reduced[i]) / (r ** p - 1))
        if err > max_error:
            max_error = err
    return max_error

# —————————————————————— Метод Рунге–Кутты ——————————————————————
class RungeKutta:
    def __init__(self, x0, y0, z0, finish):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.finish = finish

    def Get(self, h, dz, dy):
        xk, yk, zk = self.x0, self.y0, self.z0
        self.Xk = [xk]
        self.Yk = [yk]
        self.Zk = [zk]

        while abs(xk - self.finish) > abs(h) / 10:
            L1 = h * dz(xk, yk, zk)
            K1 = h * dy(xk, yk, zk)

            L2 = h * dz(xk + h / 2, yk + K1 / 2, zk + L1 / 2)
            K2 = h * dy(xk + h / 2, yk + K1 / 2, zk + L1 / 2)

            L3 = h * dz(xk + h / 2, yk + K2 / 2, zk + L2 / 2)
            K3 = h * dy(xk + h / 2, yk + K2 / 2, zk + L2 / 2)

            L4 = h * dz(xk + h, yk + K3, zk + L3)
            K4 = h * dy(xk + h, yk + K3, zk + L3)

            yk += (K1 + 2 * K2 + 2 * K3 + K4) / 6
            zk += (L1 + 2 * L2 + 2 * L3 + L4) / 6
            xk += h

            self.Xk.append(xk)
            self.Yk.append(yk)
            self.Zk.append(zk)

        return self.Yk

# —————————————————————— Метод стрельбы ——————————————————————
class Shooting:
    def __init__(self, mu1, mu2):
        self.mu1 = mu1
        self.mu2 = mu2
        self.all_paths = []  # Добавлено

    def Get(self, h):
        a, b = 2, 1
        z_b = -1
        z_a = lambda mu: (2 * mu - 4) / 4
        eps = 1e-6

        dy = lambda x, y, z: z
        dz = lambda x, y, z: 2 * y / (x**2 * (x + 1))

        mu1 = self.mu1
        mu2 = self.mu2

        rk1 = RungeKutta(a, mu1, z_a(mu1), b)
        rk2 = RungeKutta(a, mu2, z_a(mu2), b)

        rk1.Get(h, dz, dy)
        rk2.Get(h, dz, dy)

        self.all_paths.append((rk1.Xk, rk1.Yk))  # Добавлено
        self.all_paths.append((rk2.Xk, rk2.Yk))  # Добавлено

        phi1 = rk1.Zk[-1] - z_b
        phi2 = rk2.Zk[-1] - z_b

        while abs(phi2) > eps:
            dphi = (phi2 - phi1) / (mu2 - mu1)
            mu1, phi1 = mu2, phi2
            mu2 = mu2 - phi2 / dphi
            rk = RungeKutta(a, mu2, z_a(mu2), b)
            rk.Get(h, dz, dy)
            phi2 = rk.Zk[-1] - z_b

            self.all_paths.append((rk.Xk, rk.Yk))  # Добавлено

        self.Xi = rk.Xk
        self.Yi = rk.Yk
        self.mu = mu2
        return self.Yi


# —————————————————————— Метод конечных разностей ——————————————————————
class TridiagonalSolver:
    def __init__(self, A, d):
        n = len(d)
        p = [0] * n
        q = [0] * n
        x = [0] * n

        p[0] = -A[0][1] / A[0][0]
        q[0] = d[0] / A[0][0]

        for i in range(1, n - 1):
            denom = A[i][i] + A[i][i - 1] * p[i - 1]
            p[i] = -A[i][i + 1] / denom
            q[i] = (d[i] - A[i][i - 1] * q[i - 1]) / denom

        x[n - 1] = (d[n - 1] - A[n - 1][n - 2] * q[n - 2]) / (A[n - 1][n - 1] + A[n - 1][n - 2] * p[n - 2])

        for i in range(n - 2, -1, -1):
            x[i] = p[i] * x[i + 1] + q[i]

        self.x = x

    def Ans(self):
        return self.x

class FiniteDifferenceMethod:
    def Get(self, n):
        x0, xn = 1, 2
        h = (xn - x0) / n
        X = [x0 + i * h for i in range(n + 1)]

        A = [[0] * (n + 1) for _ in range(n + 1)]
        b = [0] * (n + 1)

        A[0][0] = 1
        b[0] = 2

        for i in range(1, n):
            x = X[i]
            A[i][i - 1] = 1 / h**2
            A[i][i] = -2 / h**2 + (-2 / ((x + 1) * x**2))
            A[i][i + 1] = 1 / h**2
            b[i] = 0

        A[n][n] = 1
        b[n] = 1.5

        solver = TridiagonalSolver(A, b)
        self.Xk = X
        return solver.Ans()

# —————————————————————— Основная часть ——————————————————————
def f_exact(x):
    return 1 / x + 1

# Метод стрельбы
sh = Shooting(mu1=1, mu2=3)
Y_sh = sh.Get(-0.05)

X_sh = sh.Xi

X = [1 + i * 0.01 for i in range(101)]
Y = [f_exact(x) for x in X]

PlotShootingIterations(sh.all_paths, X, Y)

print("Метод стрельбы:")
print(f" y(2) = {f_exact(2)}    y(2) ≈ mu = {sh.mu}")
print(f"Ошибка Рунге–Ромберга: {RungeRombergMet(sh, -0.05, 2, 2)}")
Plot2(X_sh, Y_sh, X, Y, "Метод стрельбы")

# Метод конечных разностей
fd = FiniteDifferenceMethod()
Y_fd = fd.Get(40)
X_fd = fd.Xk

print("\nМетод конечных разностей:")
print(f"Ошибка Рунге–Ромберга: {RungeRombergMet(fd, 40, 2, 1)}")
Plot2(X_fd, Y_fd, X, Y, "Метод конечных разностей")
