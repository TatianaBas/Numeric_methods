"""
Ввод:

3

5 5 3
5 -4 1
3 1 2

0.0001

Вывод:

Eigenvalues: [8.7055, -6.2394, 0.5339]
Eigenvectors:
[0.8303, -0.4152, -0.3717]
[0.3602, 0.9088, -0.2105]
[0.4252, 0.0409, 0.9042]
Iterations: 5

Reconstructed matrix from eigenvectors and eigenvalues:
  5.0000   5.0000   3.0000
  5.0000  -4.0000   1.0000
  3.0000   1.0000   2.0000
  
"""

import math

def L2_norm(X):
    """Find L2-norm"""
    norm = 0
    for i in range(len(X)):
        for j in range(i + 1, len(X)):
            norm += X[i][j] * X[i][j]
    return math.sqrt(norm)

def find_max_upper_element(X):
    n = len(X)
    i_max, j_max = 0, 1
    max_elem = abs(X[0][1])
    for i in range(n):
        for j in range(i + 1, n):
            if abs(X[i][j]) > max_elem:
                max_elem = abs(X[i][j])
                i_max = i
                j_max = j
    return i_max, j_max

def is_symmetric(A):
    """Check if matrix A is symmetric"""
    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            if A[i][j] != A[j][i]:
                return False
    return True

def rotation_method(A, eps):
    """Find eigenvalue and eigenvector using rotation method"""

    if not is_symmetric(A):
        raise ValueError("Matrix A must be symmetric")
    n = len(A)
    A_i = [row[:] for row in A]
    eigen_vectors = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    iterations = 0

    while L2_norm(A_i) > eps:
        i_max, j_max = find_max_upper_element(A_i)
        if A_i[i_max][i_max] - A_i[j_max][j_max] == 0:
            phi = math.pi / 4
        else:
            phi = 0.5 * math.atan(2 * A_i[i_max][j_max] / (A_i[i_max][i_max] - A_i[j_max][j_max]))

        cos_phi = math.cos(phi)
        sin_phi = math.sin(phi)

        for k in range(n):
            if k != i_max and k != j_max:
                A_ik = A_i[i_max][k] * cos_phi + A_i[j_max][k] * sin_phi
                A_jk = -A_i[i_max][k] * sin_phi + A_i[j_max][k] * cos_phi
                A_i[i_max][k] = A_ik
                A_i[j_max][k] = A_jk
                A_i[k][i_max] = A_ik  
                A_i[k][j_max] = A_jk

        A_ii = A_i[i_max][i_max] * cos_phi**2 + 2 * A_i[i_max][j_max] * cos_phi * sin_phi + A_i[j_max][j_max] * sin_phi**2
        A_jj = A_i[i_max][i_max] * sin_phi**2 - 2 * A_i[i_max][j_max] * cos_phi * sin_phi + A_i[j_max][j_max] * cos_phi**2
        A_ij = 0  

        A_i[i_max][i_max] = A_ii
        A_i[j_max][j_max] = A_jj
        A_i[i_max][j_max] = A_ij
        A_i[j_max][i_max] = A_ij

        for k in range(n):
            eig_ik = eigen_vectors[k][i_max] * cos_phi + eigen_vectors[k][j_max] * sin_phi
            eig_jk = -eigen_vectors[k][i_max] * sin_phi + eigen_vectors[k][j_max] * cos_phi
            eigen_vectors[k][i_max] = eig_ik
            eigen_vectors[k][j_max] = eig_jk

        iterations += 1

    eigen_values = [A_i[i][i] for i in range(n)]
    return eigen_values, eigen_vectors, iterations


def create_matrix(eigen_values, eigen_vectors):
    """Check: create matrix A from eigenvalues and eigenvectors"""
    n = len(eigen_values)
    A = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                A[i][j] += eigen_vectors[i][k] * eigen_values[k] * eigen_vectors[j][k]
    return A

def print_matrix(matrix):
    for row in matrix:
        print(' '.join(f'{val:8.4f}' for val in row))


if __name__ == '__main__':
    n = int(input('Enter the size of matrix: '))
    print('Enter matrix A')
    A = [list(map(float, input().split())) for _ in range(n)]
    eps = float(input('Enter epsilon: '))
    decimal_places = len(str(eps).split('.')[-1]) if '.' in str(eps) else 0

    eig_values, eig_vectors, iters = rotation_method(A, eps)
    eig_values = [round(x, decimal_places) for x in eig_values]
    print('Eigenvalues:', eig_values)
    print('Eigenvectors:')
    for row in eig_vectors:
        row = [round(x, decimal_places) for x in row]
        print(row)
    print('Iterations:', iters)

    reconstructed_matrix = create_matrix(eig_values, eig_vectors)
    print("\nReconstructed matrix from eigenvectors and eigenvalues:")
    print_matrix(reconstructed_matrix)