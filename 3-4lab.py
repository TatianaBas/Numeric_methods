"""
Точка совпадает с правой границей отрезка
Значение левосторонней производной: 0.5753600000000001
Значение правосторонней производной: 0.44628
Значение производной второго порядка точности: 0.51082
Значение второй производной: -0.25816000000000017
"""

def get_i(x, x0):
    for i in range(len(x)):
        if x[i] == x0:
            return i
    return -1  


def first_der(x, y, index_x, x0):
    if x0 in x and x0 != x[0] and x0 != x[-1]:
        print("\nТочка совпадает с правой границей отрезка")
        dxLeft = (y[index_x] - y[index_x - 1]) / (x[index_x] - x[index_x - 1])
        print("Значение левосторонней производной:", dxLeft)
        dxRight = (y[index_x + 1] - y[index_x]) / (x[index_x + 1] - x[index_x])
        print("Значение правосторонней производной:", dxRight)
        dxSecondAccuracy = (dxLeft + dxRight) / 2
        print("Значение производной второго порядка точности:", dxSecondAccuracy)
    else:
        print("\nТочка находится внутри отрезка. Найдем производную с первым порядком точности:")
        dx = (y[index_x] - y[index_x - 1]) / (x[index_x] - x[index_x - 1])
        print("Значение производной:", dx)

        if index_x < len(x) - 1:
            dxLeft = (y[index_x] - y[index_x - 1]) / (x[index_x] - x[index_x - 1])
            dxRight = (y[index_x + 1] - y[index_x]) / (x[index_x + 1] - x[index_x])
            dxSecondAccuracy = dxLeft + (dxRight - dxLeft) / (x[index_x + 1] - x[index_x - 1]) * (2 * x[index_x] - x[index_x - 1] - x[index_x])
            print("Значение производной второго порядка точности:", dxSecondAccuracy)


def second_der(x, y, index_x):
    if index_x < len(x) - 1:
        dxdx = (2 * ((y[index_x + 1] - y[index_x]) / (x[index_x + 1] - x[index_x]) -
                     (y[index_x] - y[index_x - 1]) / (x[index_x] - x[index_x - 1])) /
                (x[index_x + 1] - x[index_x - 1]))
        print("Значение второй производной:", dxdx)
    else:
        print("Невозможно вычислить производную второго порядка")


# Входные данные
x = [1.0, 1.5, 2.0, 2.5, 3.0]
y = [0.0, 0.40547, 0.69315, 0.91629, 1.0986]
x0 = 2.0

index_x = get_i(x, x0)
if index_x != -1:
    first_der(x, y, index_x, x0)
    second_der(x, y, index_x)
else:
    print("Error")
