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

if __name__ == '__main__':
    n = int(input('Enter the size of matrix: '))
    print('Enter matrix A')
    A = [list(map(float, input().split())) for _ in range(n)]
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
    print_matrix(inverse_matrix(LU))

