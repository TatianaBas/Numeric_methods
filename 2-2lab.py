"""
Метод Ньютона
Введите x1: 6
Введите x2: 1
Введите точность: 0.001
Решение системы методом Ньютона:
x = 5.929, y = 1.251
Метод простых итераций
Введите x1: 6
Введите x2: 1
Введите точность: 0.001
Решение системы методом простых итераций:
x = 5.929, y = 1.251
"""

import matplotlib.pyplot as plt
import matplotlib
import math
from decimal import Decimal

matplotlib.use('TkAgg')


def get_decimal_places(num):
    d = Decimal(str(num)) 
    return abs(d.as_tuple().exponent)

def f1(x, a):
    return (x[0] ** 2 + a ** 2) * x[1] - a ** 3


def f2(x, a):
    return (x[0] - a / 2) ** 2 + (x[1] - a / 2) ** 2 - a ** 2


def F(x, a):
    return [f1(x, a), f2(x, a)]


def jacobian(x, a):
    df1_dx0 = 2 * x[0] * x[1]
    df1_dx1 = x[0] ** 2 + a ** 2
    df2_dx0 = 2 * (x[0] - a / 2)
    df2_dx1 = 2 * (x[1] - a / 2)
    return [[df1_dx0, df1_dx1], [df2_dx0, df2_dx1]]


def vector_norm(v):
    return math.sqrt(sum(x ** 2 for x in v))


def solve_linear_2x2(A, b):
    a11, a12 = A[0]
    a21, a22 = A[1]
    b1, b2 = b

    det = a11 * a22 - a12 * a21
    if abs(det) < 1e-12:
        raise Exception("Матрица Якоби вырождена")

    x = (b1 * a22 - b2 * a12) / det
    y = (a11 * b2 - a21 * b1) / det
    return [x, y]


def newton_method(x0, a, tol, max_iter=1000):
    x = x0[:]
    for _ in range(max_iter):
        J = jacobian(x, a)
        F_val = F(x, a)
        delta_x = solve_linear_2x2(J, [-f for f in F_val])
        x_new = [x[i] + delta_x[i] for i in range(2)]

        if vector_norm(delta_x) < tol:
            return x_new

        # Ограничение на положительность координат
        x_new[0] = math.copysign(abs(x_new[0]), x0[0])
        x_new[1] = math.copysign(abs(x_new[1]), x0[1])

        x = x_new
    raise Exception("Метод Ньютона не сошелся")


def g1(x, a):
    val = a ** 2 - (x[1] - a / 2) ** 2
    if val < 0:
        raise Exception("Недопустимое значение под корнем в g1")
    return math.copysign(math.sqrt(val), x[0]) + a / 2

def g2(x, a):
    denom = x[0] ** 2 + a ** 2
    if denom == 0:
        raise Exception("Деление на ноль в g2")
    return a ** 3 / denom

def g1_der1(x, a):
    val = a ** 2 - (x[1] - a / 2) ** 2
    if val <= 0:
        raise Exception("Недопустимое значение под корнем в g1'")
    return 0  # g1 зависит только от x[1]

def g1_der2(x, a):
    val = a ** 2 - (x[1] - a / 2) ** 2
    if val <= 0:
        raise Exception("Недопустимое значение под корнем в g1'")
    return -(x[1] - a / 2) / math.sqrt(val)

def g2_der1(x, a):
    denom = (x[0] ** 2 + a ** 2) ** 2
    if denom == 0:
        raise Exception("Деление на ноль в g2'")
    return -2 * a ** 3 * x[0] / denom

def g2_der2(x, a):
    return 0  # g2 не зависит от x[1]

def G(x, a):
    return [g1(x, a), g2(x, a)]

def vector_norm(v):
    return math.sqrt(sum(val ** 2 for val in v))

def convergence_condition(x, a):
    a1 = abs(g1_der1(x, a)) + abs(g1_der2(x, a))
    a2 = abs(g2_der1(x, a)) + abs(g2_der2(x, a))
    q = max(a1, a2)
    return q < 1

def simple_iteration_method(x0, a, tol, max_iter=1000):
    x = x0[:]
    for i in range(max_iter):
        x_new = G(x, a)
        diff = [x_new[i] - x[i] for i in range(len(x))]
        if vector_norm(diff) < tol:
            return x_new
        x = x_new
    raise Exception("Метод простой итерации не сошелся")



def draw_graph(a):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_title("Система нелинейных уравнений")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)
    ax.set_xlim(-3, 7)
    ax.set_ylim(-3, 7)

    x_vals = [i * 0.1 for i in range(-30, 70)]
    y_vals = [i * 0.1 for i in range(-30, 70)]

    Z1 = [[f1([x, y], a) for x in x_vals] for y in y_vals]
    Z2 = [[f2([x, y], a) for x in x_vals] for y in y_vals]

    # Отрисовка уровней без меток
    ax.contour(x_vals, y_vals, Z1, levels=[0], colors='blue', linewidths=2)
    ax.contour(x_vals, y_vals, Z2, levels=[0], colors='red', linewidths=2)

    # Ручное добавление легенды
    legend_elements = [
        Line2D([0], [0], color='blue', lw=2, label='f1(x, y) = 0'),
        Line2D([0], [0], color='red', lw=2, label='f2(x, y) = 0')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    plt.show()



# Основной блок
a = 4
draw_graph(a)

print('Метод Ньютона')

x1 = float(input("Введите x1: "))
x2 = float(input("Введите x2: "))
x0 = [x1, x2]
tol = float(input('Введите точность: '))

solution = newton_method(x0, a, tol)
dec = get_decimal_places(tol)
print("Решение системы методом Ньютона:")
print(f"x = {round(solution[0], dec)}, y = {round(solution[1], dec)}")

print('Метод простых итераций')

x1 = float(input("Введите x1: "))
x2 = float(input("Введите x2: "))

if not convergence_condition(x0, a):
    print("Условие сходимости не выполнено — метод может не сойтись.")

x0 = [x1, x2]
tol = float(input('Введите точность: '))
dec = get_decimal_places(tol)
solution_simple_iteration = simple_iteration_method(x0, a, tol)
print("Решение системы методом простых итераций:")
print(f"x = {round(solution_simple_iteration[0], dec)}, y = {round(solution_simple_iteration[1], dec)}")