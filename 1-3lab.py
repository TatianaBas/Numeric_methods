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
[1.0000054874029314, 2.000002710606859, 2.999996803157834, -1.9999958525730877]
Iterations: 10

Seidel method
[0.9999824807788378, 1.9999990888927066, 2.99999678635388, -2.0000005516184944]
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

def solve_iterative(A, b, eps):
    """
    Метод простых итераций для решения Ax = b.
    """
    n = len(A)
    alpha = [[0] * n for _ in range(n)]
    beta = [0] * n

    # Преобразуем уравнение
    for i in range(n):
        for j in range(n):
            if i == j:
                alpha[i][j] = 0
            else:
                alpha[i][j] = -A[i][j] / A[i][i]
        beta[i] = b[i] / A[i][i]

    # Итерации
    iterations = 0
    cur_x = beta[:]
    converge = False
    while not converge:
        prev_x = cur_x[:]
        cur_x = [sum(alpha[i][j] * prev_x[j] for j in range(n)) + beta[i] for i in range(n)]
        iterations += 1
        converge = L2_norm([prev_x[i] - cur_x[i] for i in range(n)]) <= eps

    return cur_x, iterations

def seidel_multiplication(alpha, x, beta):
    """
    Умножение alpha * x + beta для метода Зейделя.
    """
    res = x[:]
    for i in range(len(alpha)):
        res[i] = beta[i]
        for j in range(len(alpha[i])):
            res[i] += alpha[i][j] * (res[j] if j < i else x[j])
    return res

def solve_seidel(A, b, eps):
    """
    Метод Зейделя для решения Ax = b.
    """
    n = len(A)
    alpha = [[0] * n for _ in range(n)]
    beta = [0] * n

    for i in range(n):
        for j in range(n):
            if i == j:
                alpha[i][j] = 0
            else:
                alpha[i][j] = -A[i][j] / A[i][i]
        beta[i] = b[i] / A[i][i]

    iterations = 0
    cur_x = beta[:]
    converge = False
    while not converge:
        prev_x = cur_x[:]
        cur_x = seidel_multiplication(alpha, cur_x, beta)
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

    print('Iteration method')
    x_iter, i_iter = solve_iterative(A, b, eps)
    print(x_iter)
    print('Iterations:', i_iter)
    print()

    print('Seidel method')
    x_seidel, i_seidel = solve_seidel(A, b, eps)
    print(x_seidel)
    print('Iterations:', i_seidel)
