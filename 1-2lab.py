"""
Ввод:
5

13 -5 0 0 0
-4 9 -5 0 0
0 -1 -12 -6 0
0 0 6 20 -5
0 0 0 4 5

-66 -47 -43 -74 14

Вывод:
Решение
[-7.0, -5.0, 6.000000000000001, -4.0, 6.0]

"""


def read_tridiagonal_matrix(n):
    """
    Считывает трёхдиагональную матрицу размером n x n, включая нулевые элементы.
    Возвращает матрицу A в виде списка списков.
    """
    A = []
    print("Введите матрицу построчно (через пробел):")
    for i in range(n):
        row = list(map(int, input().split()))
        if len(row) != n:
            raise ValueError(f"Ошибка: в строке {i + 1} должно быть {n} элементов.")
        A.append(row)
    return A


def solve(A, b):
    n = len(A)
    ok = True
    for i in range(1, n-1):
        if abs(A[i][i]) < abs(A[i][i-1]) + abs(A[i][i+1]):
            ok = False
            return ValueError("Достаточное условие не выполнено!")
    if(abs(A[0][0]) < abs(A[0][1])) or (abs(A[n-1][n-1]) < abs (A[n-1][n-2])):
        ok = False
        return ValueError("Достаточное условие не выполнено!")
    # Шаг 1
    v = [0 for i in range(n)]
    u = [0 for i in range(n)]
    v[0] = A[0][1] / -A[0][0]
    u[0] = b[0] / A[0][0]
    for i in range(1, n-1):
        v[i] = A[i][i+1] / (-A[i][i] - A[i][i-1] * v[i-1])
        u[i] = (A[i][i-1] * u[i-1] - b[i]) / (-A[i][i] - A[i][i-1] * v[i-1])
    v[n-1] = 0
    u[n-1] = (A[n-1][n-2] * u[n-2] - b[n-1]) / (-A[n-1][n-1] - A[n-1][n-2] * v[n-2])

    # Шаг 2
    x = [0 for i in range(n)]
    x[n-1] = u[n-1]
    for i in range(n-1, 0, -1):
        x[i-1] = v[i-1] * x[i] + u[i-1]
    return x

def check_solution(A, x):
    """
    Умножает матрицу A на вектор x.
    Возвращает вектор результата.
    """
    result = [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]
    return result


if __name__ == "__main__":
    n = int(input('Введите количество уравнений: '))
    # print('Введите ненулевые элементы трехдиагональной матрицы')
    A = read_tridiagonal_matrix(n)
    print('Введите вектор b')
    b = list(map(int, input().split()))

    print('Решение')
    x = solve(A, b)
    print(x)

    print("Проверка решения:")
    result = check_solution(A, x)
    print(f'A * x = {result}')