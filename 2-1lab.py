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

# Решение методом Ньютона и простых итераций решение уравнения sqrt(1 - x**2) - exp(x) + 0.1 = 0

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

def iteration_method():
    while(1):
        print('Введите концы отрезка')
        a, b = map(float, input().split())
        precision = float(input('Введите точность: '))
        decimal_places = get_decimal_places(precision)
        if(not(phi_derivative(a) < 1 and phi_derivative(b) < 1)):
            print("Производная функции на концах отрезка не меньше 1! Метод простых итераций неприменим.")
        else:
            x=(a-b)/2 
            break
    while(1):
        x_prev = x
        i = 0
        while (1):
            i += 1
            x = phi(x_prev)
            if abs(f(x) - f(x_prev)) < precision:
                break
            x_prev = x
        break
    print(f"Число итераций: {i}")
    return x, decimal_places

def initialX(x):  
    return second_derivative(x) * f(x) > 0


def convergence_condition(x):  # достаточное условие сходимости
    return abs(f(x) * second_derivative(x)) < abs(derivative(x)) ** 2

def check_conditions_newton(a, b):
    fact = True
    if convergence_condition(a) and convergence_condition(b):
        if initialX(a):
            x = a
        else:
            x = b
        print(f"Условие выполнено в обеих точках, выбираем ближайшую")
    elif convergence_condition(a):
        print(f"Условие не выполнено в точке {b}, берем {a}")
        x = a
    elif convergence_condition(b):
        print(f"Условие не выполнено в точке {a}, берем {b}")
        x = b
    else:
        print(f"Не выполняется условие сходимости в точках {a} и {b}")
        fact = False
    return fact, x

def NewtonsMethod():
    """Найти корень методом Ньютона."""
    counter = 0
    ans = 0
    while(1):
        print('Введите концы отрезка')
        a, b = map(float, input().split())
        precision = float(input('Введите точность: '))
        decimal_places = get_decimal_places(precision)
        fact, ans = check_conditions_newton(a,b)
        if(fact):
            break
    while (1):
        counter += 1
        new_ans = ans - f(ans)/derivative(ans)
        if (abs(new_ans-ans)<precision):
            ans = new_ans
            break
        ans = new_ans
    print(f"Число итераций: {counter}")
    return ans, decimal_places


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




if __name__ == "__main__":
    show()

    print('Метод Ньютона')
    x, dec = NewtonsMethod()
    func = f(x)
    print(f"x = {round(x, dec)}, f(x) = {round(func, dec)}")

    print('Метод простых итераций')
    x, dec = iteration_method()
    func = f(x)
    print(f"x = {round(x, dec)}, f(x) = {round(func, dec)}")

