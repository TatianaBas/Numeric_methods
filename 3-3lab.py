"""
Многочлен 1-й степени: f(x) = -0.19013 + 1.24621 * x
Многочлен 2-й степени: f(x) = -0.46450 + -0.12563 * x + 0.50809 * x^2

Сумма квадратов ошибок (1-й степени): 8.67902
Сумма квадратов ошибок (2-й степени): 2.35575
"""

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


class Approximation:
    def gauss(self, A, B):
        n = len(B)
        for i in range(n):
            max_elem = A[i][i]
            for j in range(i, n):
                A[i][j] /= max_elem
            B[i] /= max_elem

            for k in range(i + 1, n):
                factor = A[k][i]
                for j in range(i, n):
                    A[k][j] -= factor * A[i][j]
                B[k] -= factor * B[i]

        X = [0] * n
        for i in range(n - 1, -1, -1):
            X[i] = B[i]
            for j in range(i + 1, n):
                X[i] -= A[i][j] * X[j]
        return X

    def get_linear_coeffs(self, x, y):
        n = len(x)
        sx = sum(x)
        sx2 = sum(xi ** 2 for xi in x)
        sy = sum(y)
        sxy = sum(x[i] * y[i] for i in range(n))

        A = [[n, sx],
             [sx, sx2]]
        B = [sy, sxy]

        return self.gauss(A, B)

    def get_square_coeffs(self, x, y):
        n = len(x)
        sx = sum(x)
        sx2 = sum(xi ** 2 for xi in x)
        sx3 = sum(xi ** 3 for xi in x)
        sx4 = sum(xi ** 4 for xi in x)
        sy = sum(y)
        sxy = sum(x[i] * y[i] for i in range(n))
        sx2y = sum((x[i] ** 2) * y[i] for i in range(n))

        A = [[n, sx, sx2],
             [sx, sx2, sx3],
             [sx2, sx3, sx4]]
        B = [sy, sxy, sx2y]

        return self.gauss(A, B)


# Данные
x_data = [-0.9, 0.0, 0.9, 1.8, 2.7, 3.6]
y_data = [-0.36892, 0.0, 0.36892, 0.85408, 1.7856, 6.3138]

approximation = Approximation()

a0, a1 = approximation.get_linear_coeffs(x_data, y_data)
b0, b1, b2 = approximation.get_square_coeffs(x_data, y_data)

# Многочлены
f_lin = lambda x: a0 + a1 * x
f_sq = lambda x: b0 + b1 * x + b2 * x * x

print(f"Многочлен 1-й степени: f(x) = {a0:.5f} + {a1:.5f} * x")
print(f"Многочлен 2-й степени: f(x) = {b0:.5f} + {b1:.5f} * x + {b2:.5f} * x^2")

# Ошибки
error1 = sum((f_lin(x_data[i]) - y_data[i]) ** 2 for i in range(len(x_data)))
error2 = sum((f_sq(x_data[i]) - y_data[i]) ** 2 for i in range(len(x_data)))

print()
print(f"Сумма квадратов ошибок (1-й степени): {error1:.5f}")
print(f"Сумма квадратов ошибок (2-й степени): {error2:.5f}")

# График
x_val = [i / 100.0 for i in range(-200, 501)]
y_val1 = [f_lin(x) for x in x_val]
y_val2 = [f_sq(x) for x in x_val]

plt.figure(figsize=(10, 6))
plt.scatter(x_data, y_data, color='black', label='Исходные данные')
plt.plot(x_val, y_val1, label='Многочлен 1-й степени', color='blue')
plt.plot(x_val, y_val2, label='Многочлен 2-й степени', color='green')

plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Аппроксимация многочленами')
plt.grid(True)
plt.legend()
plt.show()
