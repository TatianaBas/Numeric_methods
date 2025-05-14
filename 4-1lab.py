"""
Метод Эйлера:
Абсолютная погрешность: 2.99487208989222
Погрешность (Рунге-Ромберг): 1.1233016448526003 

Метод Рунге-Кутты:
Абсолютная погрешность: 0.00017213798332704755
Погрешность (Рунге-Ромберг): 7.565625246286345e-07 

Метод Адамса:
Абсолютная погрешность: 0.017533284234813262
Погрешность (Рунге-Ромберг): 0.0006833886852875004 
"""

import math
import matplotlib.pyplot as plt
from typing import Callable
import matplotlib
matplotlib.use('TkAgg')

def linspace(a: float, b: float, points: int) -> list[float]:
    return [a + (b - a) / (points - 1) * i for i in range(points)]

# Правая часть уравнения: y'' = 2y + 4x^2 * exp(x^2)
def f(x: float, y: float, z: float) -> float:
    return z

def g(x: float, y: float, z: float) -> float:
    return 2 * y + 4 * x**2 * math.exp(x**2)

# Точное решение
def y_exact(x: float) -> float:
    return math.exp(x**2) + math.exp(x * math.sqrt(2)) + math.exp(-x * math.sqrt(2))

# Метод Эйлера
class ODE2_Euler:
    def __init__(self, x0, xn, h, y0, z0, f, g):
        self.X = linspace(x0, xn, int((xn - x0) / h + 1))
        self.Y = [y0] * len(self.X)
        self.Z = [z0] * len(self.X)
        for i in range(1, len(self.X)):
            self.Z[i] = self.Z[i-1] + h * g(self.X[i-1], self.Y[i-1], self.Z[i-1])
            self.Y[i] = self.Y[i-1] + h * f(self.X[i-1], self.Y[i-1], self.Z[i-1])

    def MAE(self, true: Callable):
        return sum([abs(self.Y[i] - true(self.X[i])) for i in range(len(self.X))])

# Метод Рунге-Кутты 4-го порядка
class ODE2_Runge_Kutta:
    def __init__(self, x0, xn, h, y0, z0, f, g):
        self.X = linspace(x0, xn, int((xn - x0) / h + 1))
        self.Y = [y0] * len(self.X)
        self.Z = [z0] * len(self.X)
        for i in range(len(self.X) - 1):
            xi, yi, zi = self.X[i], self.Y[i], self.Z[i]

            K1 = h * f(xi, yi, zi)
            L1 = h * g(xi, yi, zi)

            K2 = h * f(xi + h/2, yi + K1/2, zi + L1/2)
            L2 = h * g(xi + h/2, yi + K1/2, zi + L1/2)

            K3 = h * f(xi + h/2, yi + K2/2, zi + L2/2)
            L3 = h * g(xi + h/2, yi + K2/2, zi + L2/2)

            K4 = h * f(xi + h, yi + K3, zi + L3)
            L4 = h * g(xi + h, yi + K3, zi + L3)

            self.Y[i+1] = yi + (K1 + 2*K2 + 2*K3 + K4) / 6
            self.Z[i+1] = zi + (L1 + 2*L2 + 2*L3 + L4) / 6

    def MAE(self, true: Callable):
        return sum([abs(self.Y[i] - true(self.X[i])) for i in range(len(self.X))])

# Метод Адамса
class ODE2_Adams:
    def __init__(self, x0, xn, h, y0, z0, f, g):
        self.X = linspace(x0, xn, int((xn - x0) / h + 1))
        self.Y = [y0] * len(self.X)
        self.Z = [z0] * len(self.X)
        F = [0] * len(self.X)
        G = [0] * len(self.X)

        rk = ODE2_Runge_Kutta(x0, x0 + 3*h, h, y0, z0, f, g)
        self.Y[:4] = rk.Y[:4]
        self.Z[:4] = rk.Z[:4]

        for i in range(4):
            F[i] = f(self.X[i], self.Y[i], self.Z[i])
            G[i] = g(self.X[i], self.Y[i], self.Z[i])

        for i in range(3, len(self.X) - 1):
            self.Y[i+1] = self.Y[i] + h/24 * (55*F[i] - 59*F[i-1] + 37*F[i-2] - 9*F[i-3])
            self.Z[i+1] = self.Z[i] + h/24 * (55*G[i] - 59*G[i-1] + 37*G[i-2] - 9*G[i-3])
            F[i+1] = f(self.X[i+1], self.Y[i+1], self.Z[i+1])
            G[i+1] = g(self.X[i+1], self.Y[i+1], self.Z[i+1])

    def MAE(self, true: Callable):
        return sum([abs(self.Y[i] - true(self.X[i])) for i in range(len(self.X))])

# Метод Рунге-Ромберга
class RungeRomberg:
    def __init__(self, X, Y_h, Y_h2, p):
        self.X = X
        self.Y = [y2 + (y2 - y1) / (2**p - 1) for y1, y2 in zip(Y_h, Y_h2[::2])]

    def MAE(self, true: Callable):
        return sum([abs(self.Y[i] - true(self.X[i])) for i in range(len(self.X))])

# Построение графика
def plot(solvers, labels, title):
    x = solvers[0].X
    y_exact_vals = [y_exact(xi) for xi in x]

    plt.figure(figsize=(10, 6))
    plt.plot(x, y_exact_vals, label="Точное решение", linewidth=2)

    for solver, label in zip(solvers, labels):
        plt.plot(solver.X, solver.Y, '--', label=label)

    plt.title(title)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
    plt.show()

# Главная функция
def main():
    x0, xn, h = 0, 1, 0.1
    y0, z0 = 3, 0

    euler1 = ODE2_Euler(x0, xn, h, y0, z0, f, g)
    euler2 = ODE2_Euler(x0, xn, h/2, y0, z0, f, g)
    rr_euler = RungeRomberg(euler1.X, euler1.Y, euler2.Y, 2)

    runge1 = ODE2_Runge_Kutta(x0, xn, h, y0, z0, f, g)
    runge2 = ODE2_Runge_Kutta(x0, xn, h/2, y0, z0, f, g)
    rr_runge = RungeRomberg(runge1.X, runge1.Y, runge2.Y, 4)

    adams1 = ODE2_Adams(x0, xn, h, y0, z0, f, g)
    adams2 = ODE2_Adams(x0, xn, h/2, y0, z0, f, g)
    rr_adams = RungeRomberg(adams1.X, adams1.Y, adams2.Y, 4)

    print("Метод Эйлера:")
    print("Абсолютная погрешность:", euler1.MAE(y_exact))
    print("Погрешность (Рунге-Ромберг):", rr_euler.MAE(y_exact), "\n")

    print("Метод Рунге-Кутты:")
    print("Абсолютная погрешность:", runge1.MAE(y_exact))
    print("Погрешность (Рунге-Ромберг):", rr_runge.MAE(y_exact), "\n")

    print("Метод Адамса:")
    print("Абсолютная погрешность:", adams1.MAE(y_exact))
    print("Погрешность (Рунге-Ромберг):", rr_adams.MAE(y_exact), "\n")

    plot([euler1, rr_euler], ["Эйлер", "Эйлер (Рунге-Роберт)"], "Метод Эйлера")
    plot([runge1, rr_runge], ["Рунге-Кутта", "Рунге-Кутт (Рунге-Роберт)"], "Метод Рунге-Кутты")
    plot([adams1, rr_adams], ["Адамс", "Адамс (Рунге-Роберт)"], "Метод Адамса")

if __name__ == "__main__":
    main()
