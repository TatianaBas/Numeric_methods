"""
Ввод:

4

-23 -7 5 2
-7 -21 4 9
9 5 -31 -8
0 1 -2 10

-26 -55 -58 -24

0.0001

Вывод:

Iteration method
[1.0, 2.0, 3.0, -2.0]
Iterations: 10

Seidel method
[1.0, 2.0, 3.0, -2.0]
Iterations: 6

"""

import math

def L2_norm(X):
    """
    Вычисляет норму вектора L2.
    """
    l2_norm = sum(x * x for x in X)
    return math.sqrt(l2_norm)

def matrix_vector_multiplication(A, x):
    """
    Умножает матрицу A на вектор x.
    """
    result = [0] * len(A)
    for i in range(len(A)):
        for j in range(len(A[i])):
            result[i] += A[i][j] * x[j]
    return result

def is_diagonally_dominant(A):
    """
    Проверяет достаточное условие единственности решения: диагональное преобладание.
    """
    n = len(A)
    for i in range(n):
        if abs(A[i][i]) <= sum(abs(A[i][j]) for j in range(n) if i != j):
            return False
    return True

def is_diagonally_dominant2(A):
    """
    Проверяет достаточное условие единственности решения: диагональное преобладание по столбцам.
    """
    n = len(A)
    for i in range(n):
        if abs(A[i][i]) <= sum(abs(A[j][i]) for j in range(n) if i != j):
            return False
    return True


def solve_iterative(A, b, eps):
    """
    Метод простых итераций для решения Ax = b.
    """
    n = len(A)
    iterations = 0
    cur_x = [b[i] / A[i][i] for i in range(n)]  # Начальное приближение
    converge = False
    
    while not converge:
        prev_x = cur_x[:]
        cur_x = [
            sum(-A[i][j] / A[i][i] * prev_x[j] for j in range(n) if i != j) + b[i] / A[i][i]
            for i in range(n)
        ]
        iterations += 1
        converge = L2_norm([prev_x[i] - cur_x[i] for i in range(n)]) <= eps

    return cur_x, iterations


def solve_seidel(A, b, eps):
    """
    Метод Зейделя для решения Ax = b без хранения alpha и beta.
    """
    n = len(A)
    iterations = 0
    cur_x = [b[i] / A[i][i] for i in range(n)]  # Начальное приближение
    converge = False
    
    while not converge:
        prev_x = cur_x[:]
        for i in range(n):
            sum_val = sum(-A[i][j] / A[i][i] * cur_x[j] for j in range(i)) + \
                      sum(-A[i][j] / A[i][i] * prev_x[j] for j in range(i + 1, n))
            cur_x[i] = sum_val + b[i] / A[i][i]
        iterations += 1
        converge = L2_norm([prev_x[i] - cur_x[i] for i in range(n)]) <= eps

    return cur_x, iterations


if __name__ == '__main__':
    n = int(input('Enter the size of matrix: '))

    print('Enter matrix A')
    A = [list(map(float, input().split())) for _ in range(n)]

    print('Enter vector b')
    b = list(map(float, input().split()))

    eps = float(input('Enter epsilon: '))

    if not is_diagonally_dominant(A) and not is_diagonally_dominant2(A):
        raise ValueError("Матрица не удовлетворяет достаточному условию единственности решения (диагональное преобладание).")
    decimal_places = len(str(eps).split('.')[-1]) if '.' in str(eps) else 0
    print('Iteration method')
    x_iter, i_iter = solve_iterative(A, b, eps)
    x_iter = [round(x, decimal_places) for x in x_iter]
    print(x_iter)
    print('Iterations:', i_iter)
    print()

    print('Seidel method')
    x_seidel, i_seidel = solve_seidel(A, b, eps)
    x_seidel = [round(x, decimal_places) for x in x_seidel]
    print(x_seidel)
    print('Iterations:', i_seidel)
