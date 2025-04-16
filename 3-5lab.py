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
    x = frange(x1, xk + h / 2, h)
    result = 0
    for i in x:
        if i <= x1 or i >= xk:
            result += f(i) / 2
        else:
            result += f(i)
    return result * h


def simpson(f, x1, xk, h):
    N = int((xk - x1) / h)
    if N % 2 != 0:
        N += 1  # Симпсон требует чётное число отрезков
        h = (xk - x1) / N
    result = f(x1) + f(xk)
    for i in range(1, N):
        xi = x1 + i * h
        result += 4 * f(xi) if i % 2 != 0 else 2 * f(xi)
    return result * h / 3


# Исходные данные
y = lambda x: x / (3 * x + 4) ** 3
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
