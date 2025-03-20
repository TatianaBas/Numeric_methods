"""
Ввод:

3

5 -5 -6
-1 -8 -5
2 7 -3

0.0001

Вывод:

Eigenvalues: [(-5.415796256962388+5.83122959589451j), (-5.415796256962388-5.83122959589451j), 4.831521398091571]
  
"""


import math

def sign(x):
    return -1 if x < 0 else 1 if x > 0 else 0

def L2_norm(vec):
    ans = sum(num * num for num in vec)
    return math.sqrt(ans)

def get_householder_matrix(A, col_num):
    n = len(A)
    v = [0] * n
    a = [A[i][col_num] for i in range(n)]
    v[col_num] = a[col_num] + sign(a[col_num]) * L2_norm(a[col_num:])
    for i in range(col_num + 1, n):
        v[i] = a[i]
    
    norm_v = sum(v[i] ** 2 for i in range(n))
    if norm_v == 0:
        return [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    
    H = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            H[i][j] -= 2 * v[i] * v[j] / norm_v
    
    return H

def matrix_mult(A, B):
    n = len(A)
    C = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            C[i][j] = sum(A[i][k] * B[k][j] for k in range(n))
    return C

def QR_decomposition(A):
    n = len(A)
    Q = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    A_i = [row[:] for row in A]
    
    for i in range(n - 1):
        H = get_householder_matrix(A_i, i)
        Q = matrix_mult(Q, H)
        A_i = matrix_mult(H, A_i)
    
    return Q, A_i

def get_roots(A, i):
    a11, a12 = A[i][i], A[i][i + 1] if i + 1 < len(A) else 0
    a21, a22 = A[i + 1][i] if i + 1 < len(A) else 0, A[i + 1][i + 1] if i + 1 < len(A) else 0
    D = (a11 + a22) ** 2 - 4 * (a11 * a22 - a12 * a21)
    if D < 0:
        real_part = (a11 + a22) / 2
        imag_part = math.sqrt(-D) / 2
        return (real_part + imag_part * 1j, real_part - imag_part * 1j)
    else:
        return ((a11 + a22 + math.sqrt(D)) / 2, (a11 + a22 - math.sqrt(D)) / 2)

def is_complex(A, i, eps):
    Q, R = QR_decomposition(A)
    A_next = matrix_mult(R, Q)
    lambda1 = get_roots(A, i)
    lambda2 = get_roots(A_next, i)
    return abs(lambda1[0] - lambda2[0]) <= eps and abs(lambda1[1] - lambda2[1]) <= eps

def get_eigen_value(A, i, eps):
    A_i = [row[:] for row in A]
    while True:
        Q, R = QR_decomposition(A_i)
        A_i = matrix_mult(R, Q)
        if L2_norm([A_i[k][i] for k in range(i + 1, len(A))]) <= eps:
            return A_i[i][i], A_i
        elif L2_norm([A_i[k][i] for k in range(i + 2, len(A))]) <= eps and is_complex(A_i, i, eps):
            return get_roots(A_i, i), A_i

def get_eigen_values_QR(A, eps):
    n = len(A)
    A_i = [row[:] for row in A]
    eigen_values = []
    i = 0
    while i < n:
        cur_eigen_values, A_i_plus_1 = get_eigen_value(A_i, i, eps)
        if isinstance(cur_eigen_values, tuple):
            eigen_values.extend(cur_eigen_values)
            i += 2
        else:
            eigen_values.append(cur_eigen_values)
            i += 1
        A_i = A_i_plus_1
    return eigen_values

if __name__ == '__main__':
    n = int(input('Enter the size of matrix: '))
    print('Enter matrix A')
    A = [list(map(float, input().split())) for _ in range(n)]
    eps = float(input('Enter epsilon: '))
    eig_values = get_eigen_values_QR(A, eps)
    print('Eigenvalues:', eig_values)
