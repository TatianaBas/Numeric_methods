"""
Lagrange interpolation
Points A
Polynom
L(x) = + -0.00*(x-0.39)(x-0.79)(x-1.18) + 3.42*(x-0.00)(x-0.79)(x-1.18) + -8.26*(x-0.00)(x-0.39)(x-1.18) + 6.64*(x-0.00)(x-0.39)(x-0.79)
Abs error for test point = 0.023571856732751417
Points B
Polynom
L(x) = + -0.00*(x-0.39)(x-1.05)(x-1.18) + 2.05*(x-0.00)(x-1.05)(x-1.18) + -19.31*(x-0.00)(x-0.39)(x-1.18) + 19.93*(x-0.00)(x-0.39)(x-1.05)
Abs error for test point = 0.08292780909109965

Newton interpolation
Points A
Polynom
P(x) = 0.00 + (x-0.00)*1.05 + (x-0.00)(x-0.39)*0.56 + (x-0.00)(x-0.39)(x-0.79)*1.81
Abs error for test point = 0.023571856732751306
Points B
Polynom
P(x) = 0.00 + (x-0.00)*1.05 + (x-0.00)(x-0.39)*0.92 + (x-0.00)(x-0.39)(x-1.05)*2.68
Abs error for test point = 0.08292780909109954
"""

import math

def f(x):
    return math.tan(x)

def lagrange_interpolation(x, y, test_point):
    assert len(x) == len(y)
    polynom_str = 'L(x) ='
    polynom_test_value = 0
    for i in range(len(x)):
        cur_enum_str = ''
        cur_enum_test = 1
        cur_denom = 1
        for j in range(len(x)):
            if i == j:
                continue
            cur_enum_str += f'(x-{x[j]:.2f})'
            cur_enum_test *= (test_point[0] - x[j])
            cur_denom *= (x[i] - x[j])

        coef = y[i] / cur_denom

        if i == 0:
            polynom_str += f' {coef:.2f}*' + cur_enum_str
        else:
            sign = ' + ' if coef >= 0 else ' - '
            polynom_str += sign + f'{abs(coef):.2f}*' + cur_enum_str

        polynom_test_value += y[i] * cur_enum_test / cur_denom

    abs_error = abs(polynom_test_value - test_point[1])
    return polynom_str, polynom_test_value, abs_error


def newton_interpolation(x, y, test_point):
    assert len(x) == len(y)
    n = len(x)
    coefs = [y[i] for i in range(n)]
    for i in range(1, n):
        for j in range(n - 1, i - 1, -1):
            coefs[j] = float(coefs[j] - coefs[j - 1]) / float(x[j] - x[j - i])

    polynom_str = 'P(x) = '
    polynom_test_value = 0

    cur_multipliers_str = ''
    cur_multipliers = 1
    for i in range(n):
        polynom_test_value += cur_multipliers * coefs[i]
        if i == 0:
            polynom_str += f'{coefs[i]:.2f}'
        else:
            polynom_str += ' + ' + cur_multipliers_str + '*' + f'{coefs[i]:.2f}'

        cur_multipliers *= (test_point[0] - x[i])
        cur_multipliers_str += f'(x-{x[i]:.2f})'
    return polynom_str, polynom_test_value, abs(polynom_test_value - test_point[1])

if __name__ == '__main__':

    x_a = [0, math.pi / 8, 2 * math.pi / 8, 3 * math.pi / 8]  # a)
    x_b = [0, math.pi / 8, math.pi / 3, 3 * math.pi / 8]      # b)

    y_a = [f(x) for x in x_a]
    y_b = [f(x) for x in x_b]

    x_test = 3 * math.pi / 16
    y_test = f(x_test)

    print('Lagrange interpolation')
    print('Points A')
    lagrange_polynom_a, lagrange_value_a, lagrange_error_a = lagrange_interpolation(x_a, y_a, (x_test, y_test))
    print('Polynom')
    print(lagrange_polynom_a)
    print('Interpolated value=', lagrange_value_a)
    print('Abs error for test point =', lagrange_error_a)

    print('Points B')
    lagrange_polynom_b, lagrange_value_b, lagrange_error_b = lagrange_interpolation(x_b, y_b, (x_test, y_test))
    print('Polynom')
    print(lagrange_polynom_b)
    print('Interpolated value=', lagrange_value_b)
    print('Abs error for test point =', lagrange_error_b)
    print()

    print('Newton interpolation')
    print('Points A')
    newton_polynom_a, newton_value_a, newton_error_a = newton_interpolation(x_a, y_a, (x_test, y_test))
    print('Polynom')
    print(newton_polynom_a)
    print('Interpolated value=', newton_value_a)
    print('Abs error for test point =', newton_error_a)

    print('Points B')
    newton_polynom_b, newton_value_b, newton_error_b = newton_interpolation(x_b, y_b, (x_test, y_test))
    print('Polynom')
    print(newton_polynom_b)
    print('Interpolated value=', newton_value_b)
    print('Abs error for test point =', newton_error_b)
