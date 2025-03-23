"""
Ввод:

0 0.5

0.0001

Вывод:

Метод Ньютона
Число итераций: 5
x = 0.0915, f(x) = -0.0

Метод простых итераций
Число итераций: 5
x = 0.0915, f(x) = 0.0
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from decimal import Decimal

def get_decimal_places(num):
    d = Decimal(str(num)) 
    return abs(d.as_tuple().exponent)

def f(x):
    return math.sqrt(1 - x**2) - math.exp(x) + 0.1

def derivative(x):
    return -x/math.sqrt(1 - x**2) - math.exp(x)

def second_derivative(x):
    return -((1 + x**2) / (1 - x**2)**(3/2)) - math.exp(x)

def phi(x):
    return math.log(math.sqrt(1-x**2)+0.1)

def phi_derivative(x):
    return -x / (math.sqrt(1 - x**2) * (math.sqrt(1 - x**2) + 0.1))

def iteration_method(x, eps):
    x_prev = x
    i = 0
    while True:
        i += 1
        x = phi(x_prev)
        if abs(f(x) - f(x_prev)) < eps:
            break
        x_prev = x
    print(f"Число итераций: {i}")
    return x

def NewtonsMethod(x, eps):
    if abs(x) >= 1:
        raise ValueError("x должен быть в пределах (-1,1), так как sqrt(1 - x^2) определена только в этом интервале.")
    if(not(f(x) * second_derivative(x) < 2 * (derivative(x))**2)):
        raise ValueError("Не выполняется неравенство второй производной")
    prev = x
    x = x - f(x)/derivative(x)
    i = 1
    while abs(x - prev) > eps:
        prev = x
        x = x - f(x)/derivative(x)
        i += 1
    print(f"Число итераций: {i}")
    return x

# Строим график
def show():
    x_vals = [i / 100 for i in range(-99, 100)]
    y_vals = [f(x) for x in x_vals]

    plt.figure(figsize=(8, 6))
    plt.plot(x_vals, y_vals, label=r'$y = \sqrt{1 - x^2} - e^x + 0.1$', color='b')
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("График функции")
    plt.legend()
    plt.grid(True)

    plt.show()
    return


"""
корень: [a, b] = [0; 0.5]
"""


if __name__ == "__main__":
    show()
    print('Введите концы отрезка')
    l, r = map(float, input().split())
    eps = float(input('Введите точность: '))
    decimal_places = get_decimal_places(eps)

    print('Метод Ньютона')
    x = NewtonsMethod((l-r)*0.5, eps)
    func = f(x)
    print(f"x = {round(x, decimal_places)}, f(x) = {round(func, decimal_places)}\n")

    if(not(phi_derivative(l) < 1 and phi_derivative(r) < 1)):
        raise ValueError("Производная функции на концах отрезка не меньше 1! Метод простых итераций неприменим.")
    print('Метод простых итераций')
    x = iteration_method((l-r)*0.5, eps)
    func = f(x)
    print(f"x = {round(x, decimal_places)}, f(x) = {round(func, decimal_places)}")

