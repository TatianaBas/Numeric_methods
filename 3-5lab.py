"""
Прямоугольники:
 h1=0.5: -0.07090982378616856, h2=0.25: -0.10243925004315926, уточнённый: -0.1129490587954895, ошибка: 0.010509808752330235

Трапеции:
 h1=0.5: -0.2637685074244359, h2=0.25: -0.16733916560530224, уточнённый: -0.135196051665591, ошибка: 0.032143113939711226

Симпсон:
 h1=0.5: -0.18551058521508848, h2=0.25: -0.13519605166559104, уточнённый: -0.13184174942895788, ошибка: 0.003354302236633161
"""

def frange(start, stop, step):
    values = []
    while start < stop:
        values.append(start)
        start += step
    return values


def rectangles(f, x1, xk, h):
    x = frange(x1, xk, h)
    return sum(f((i + i + h) / 2) * h for i in x)


def trapeze(f, x1, xk, h):
    N = int((xk - x1) / h)
    result = (f(x1) + f(xk)) / 2 
    for i in range(1, N):
        xi = x1 + i * h
        result += f(xi)
    return result * h


def simpson(f, x1, xk, h):
    N = int((xk - x1) / h)
    if N % 2 != 0:
        N += 1 
        h = (xk - x1) / N
    result = f(x1) + f(xk)
    for i in range(1, N):
        xi = x1 + i * h
        result += 4 * f(xi) if i % 2 != 0 else 2 * f(xi)
    return result * h / 3


# Исходные данные
y = lambda x: x / (3 * x + 4) ** 2
x0 = -1
xk = 1
h1 = 0.5
h2 = 0.25

# Прямоугольники
r1 = rectangles(y, x0, xk, h1)
r2 = rectangles(y, x0, xk, h2)
Sr = r2 + (r2 - r1) / (2**2 - 1)
print(f"Прямоугольники:\n h1={h1}: {r1}, h2={h2}: {r2}, уточнённый: {Sr}, ошибка: {abs(Sr - r2)}\n")

# Трапеции
t1 = trapeze(y, x0, xk, h1)
t2 = trapeze(y, x0, xk, h2)
St = t2 + (t2 - t1) / (2**2 - 1)
print(f"Трапеции:\n h1={h1}: {t1}, h2={h2}: {t2}, уточнённый: {St}, ошибка: {abs(St - t2)}\n")

# Симпсон
s1 = simpson(y, x0, xk, h1)
s2 = simpson(y, x0, xk, h2)
Ss = s2 + (s2 - s1) / (2**4 - 1)
print(f"Симпсон:\n h1={h1}: {s1}, h2={h2}: {s2}, уточнённый: {Ss}, ошибка: {abs(Ss - s2)}\n")
