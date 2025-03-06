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
    return math.sqrt(sum(x * x for x in X))

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
    decimal_places = len(str(eps).split('.')[-1]) if '.' in str(eps) else 0
    iterations = 0
    cur_x = [b[i] / A[i][i] for i in range(n)]  # Начальное приближение
    converge = False
    
    while not converge:
        prev_x = cur_x[:]
        cur_x = [
            round(sum(-A[i][j] / A[i][i] * prev_x[j] for j in range(n) if i != j) + b[i] / A[i][i], decimal_places)
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
    decimal_places = len(str(eps).split('.')[-1]) if '.' in str(eps) else 0
    iterations = 0
    cur_x = [b[i] / A[i][i] for i in range(n)]  # Начальное приближение
    converge = False
    
    while not converge:
        prev_x = cur_x[:]
        for i in range(n):
            sum_val = sum(-A[i][j] / A[i][i] * cur_x[j] for j in range(i)) + \
                      sum(-A[i][j] / A[i][i] * prev_x[j] for j in range(i + 1, n))
            cur_x[i] = round(sum_val + b[i] / A[i][i], decimal_places)
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
    
    print('Iteration method')
    x_iter, i_iter = solve_iterative(A, b, eps)
    print(x_iter)
    print('Iterations:', i_iter)
    print()
    
    print('Seidel method')
    x_seidel, i_seidel = solve_seidel(A, b, eps)
    print(x_seidel)
    print('Iterations:', i_seidel)

