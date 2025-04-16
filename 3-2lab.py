import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

class CubicSpline:
    def __init__(self, x, f):
        self.x = x
        self.f = f
        self.n = len(x)
        self.h = [x[i+1] - x[i] for i in range(self.n - 1)]
        self.a, self.b, self.c, self.d = self.compute_coeffs()

    def compute_coeffs(self):
        n = self.n
        h = self.h
        alpha = [0.0] * (n - 1)
        for i in range(1, n - 1):
            alpha[i] = (3 / h[i]) * (self.f[i+1] - self.f[i]) - (3 / h[i-1]) * (self.f[i] - self.f[i-1])

        l = [1.0] + [0.0] * (n - 1)
        mu = [0.0] * n
        z = [0.0] * n

        for i in range(1, n - 1):
            l[i] = 2 * (self.x[i+1] - self.x[i-1]) - h[i-1] * mu[i-1]
            mu[i] = h[i] / l[i]
            z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i]

        l[n-1] = 1.0
        z[n-1] = 0.0
        c = [0.0] * n
        b = [0.0] * (n - 1)
        d = [0.0] * (n - 1)
        a = self.f[:-1]

        for j in range(n - 2, -1, -1):
            c[j] = z[j] - mu[j] * c[j+1]
            b[j] = (self.f[j+1] - self.f[j]) / h[j] - h[j] * (c[j+1] + 2 * c[j]) / 3
            d[j] = (c[j+1] - c[j]) / (3 * h[j])

        return a, b, c, d

    def evaluate(self, x_val):
        for i in range(self.n - 1):
            if self.x[i] <= x_val <= self.x[i+1]:
                dx = x_val - self.x[i]
                return self.a[i] + self.b[i]*dx + self.c[i]*dx**2 + self.d[i]*dx**3
        return None

    def plot(self, X_star=None, f_star=None):
        xs = []
        ys = []
        step = 0.01
        x_min, x_max = self.x[0], self.x[-1]
        cur = x_min
        while cur <= x_max:
            xs.append(cur)
            ys.append(self.evaluate(cur))
            cur += step

        plt.plot(self.x, self.f, 'o', label='Исходные точки')
        plt.plot(xs, ys, '-', label='Кубический сплайн')
        if X_star is not None and f_star is not None:
            plt.plot(X_star, f_star, 's', label=f'f(x={X_star:.2f})={f_star:.4f}')
        plt.grid(True)
        plt.title("Кубический сплайн")
        plt.xlabel("x")
        plt.ylabel("f(x)")
        plt.legend()
        plt.show()


# Данные из таблицы
x_data = [0.0, 0.9, 1.8, 2.7, 3.6]
f_data = [0.0, 0.36892, 0.85408, 1.7856, 6.3138]
X_star = 1.5

spline = CubicSpline(x_data, f_data)
f_star = spline.evaluate(X_star)

print(f"Значение в x* = {X_star} : f(x*) = {f_star:.6f}")

spline.plot(X_star, f_star)
