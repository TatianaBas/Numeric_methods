"""
Ввод:
4

9 -5 -6 3
1 -7 1 0
3 -4 9 0
6 -1 9 8

-8 38 47 -8

Вывод:
LU decomposition
  9.00  -5.00  -6.00   3.00
  0.11  -6.44   1.67  -0.33
  0.33   0.36  10.40  -0.88
  0.67  -0.36   1.31   7.03
System solution
x: [3.9474596431116675e-16, -5.0, 3.0, -5.000000000000001]
det A = -4239.0
A^(-1)
  0.11  -0.15   0.13  -0.04
  0.01  -0.17   0.03  -0.00
 -0.03  -0.02   0.08   0.01
 -0.05   0.12  -0.19   0.14
A * A^(-1)
  1.00   0.00   0.00   0.00
  0.00   1.00   0.00  -0.00
  0.00   0.00   1.00   0.00
  0.00   0.00   0.00   1.00
check solution
A* x =  [-8.0, 38.0, 47.0, -8.000000000000007]

"""

import copy
import math

def LU_decompose(A):
    """
    A = LU, where L and U are stored in the same matrix.
    """
    n = len(A)
    for k in range(n):
        for i in range(k + 1, n):
            A[i][k] /= A[k][k]
            for j in range(k + 1, n):
                A[i][j] -= A[i][k] * A[k][j]
    return A

def solve_system(A, b):
    """Solves system of equations: LUx = b using combined LU matrix."""
    n = len(A)
    y = [0] * n
    for i in range(n):
        y[i] = b[i] - sum(A[i][j] * y[j] for j in range(i))

    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]
    return x

def determinant(A):
    """Calculate the determinant of matrix A."""
    return math.prod(A[i][i] for i in range(len(A)))

def inverse_matrix(A):
    """Calculate A^(-1) using combined LU matrix."""
    n = len(A)
    E = [[float(i == j) for j in range(n)] for i in range(n)]
    A_inv = [solve_system(A, E[i]) for i in range(n)]
    return list(map(list, zip(*A_inv)))

def print_matrix(A):
    for row in A:
        print(' '.join(f'{val:6.2f}' for val in row))

def multiply_matrix_vector(A, x):
    """Проверка."""
    n = len(A)
    result = [0] * n
    for i in range(n):
        result[i] = sum(A[i][j] * x[j] for j in range(n))
    return result

def multiply_matrices(A, B):
    """Перемножает две матрицы A и B."""
    n = len(A)
    result = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            result[i][j] = sum(A[i][k] * B[k][j] for k in range(n))
    return result

if __name__ == '__main__':
    n = int(input('Enter the size of matrix: '))
    print('Enter matrix A')
    A = [list(map(float, input().split())) for _ in range(n)]
    A_copy = [row[:] for row in A]
    print('Enter vector b')
    b = list(map(float, input().split()))

    print("LU decomposition")
    LU = LU_decompose(A)
    print_matrix(LU)

    print("System solution")
    x = solve_system(LU, b)
    print('x:', x)

    print("det A =", determinant(LU))

    print("A^(-1)")
    A_inv = inverse_matrix(LU)
    print_matrix(A_inv)

    print("A * A^(-1)")
    identity = multiply_matrices(A_copy, A_inv)
    print_matrix(identity)


    print("check solution")
    Ax = multiply_matrix_vector(A_copy, x)
    print("A* x = ", Ax)
