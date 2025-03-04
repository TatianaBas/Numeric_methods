"""
Ввод:

3

5 5 3
5 -4 1
3 1 2

0.0001

Вывод:

Eigenvalues: [8.705467377749228, -6.239373522935353, 0.5339061451861218]
Eigenvectors:
[0.8303244236416637, -0.4152103043947342, -0.37170116307131923]
[0.3602250084569686, 0.9088062601455725, -0.21049732730463597]
[0.425205003247041, 0.040885017361501484, 0.9041731695693218]
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


def rotation_method(A, eps):
    """Find eigenvalue and eigenvector using rotation method"""
    n = len(A)
    A_i = [row[:] for row in A]
    eigen_vectors = [[0.0] * n for _ in range(n)]
    for i in range(n):
        eigen_vectors[i][i] = 1.0
    iterations = 0

    while L2_norm(A_i) > eps:
        i_max, j_max = find_max_upper_element(A_i)
        if A_i[i_max][i_max] - A_i[j_max][j_max] == 0:
            phi = math.pi / 4
        else:
            phi = 0.5 * math.atan(2 * A_i[i_max][j_max] / (A_i[i_max][i_max] - A_i[j_max][j_max]))

        """rotation matrix"""
        U = [[0.0] * n for _ in range(n)]
        for i in range(n):
            U[i][i] = 1.0 

        U[i_max][j_max] = -math.sin(phi)
        U[j_max][i_max] = math.sin(phi)
        U[i_max][i_max] = math.cos(phi)
        U[j_max][j_max] = math.cos(phi)

        U_T = [[U[j][i] for j in range(n)] for i in range(n)]
        temp_matrix = [[sum(U_T[i][k] * A_i[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        A_i = [[sum(temp_matrix[i][k] * U[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        temp_eigen_vectors = [[sum(eigen_vectors[i][k] * U[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        eigen_vectors = temp_eigen_vectors

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

    eig_values, eig_vectors, iters = rotation_method(A, eps)
    print('Eigenvalues:', eig_values)
    print('Eigenvectors:')
    for row in eig_vectors:
        print(row)
    print('Iterations:', iters)

    reconstructed_matrix = create_matrix(eig_values, eig_vectors)
    print("\nReconstructed matrix from eigenvectors and eigenvalues:")
    print_matrix(reconstructed_matrix)