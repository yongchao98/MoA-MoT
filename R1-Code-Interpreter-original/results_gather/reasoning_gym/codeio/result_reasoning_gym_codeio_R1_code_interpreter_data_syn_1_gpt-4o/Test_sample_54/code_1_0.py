import math

def grule(n):
    m = (n + 1) // 2
    e1 = n * (n + 1)

    x = []
    w = []

    for i in range(m):
        t = (4 * i + 3) * math.pi / (4 * n + 2)
        x0 = (1.0 - (1.0 - 1.0 / n) / (8.0 * n * n)) * math.cos(t)
        pkm1 = 1.0
        pk = x0

        for k in range(2, n + 1):
            t1 = x0 * pk
            pkp1 = t1 - pkm1 - (t1 - pkm1) / k + t1
            pkm1 = pk
            pk = pkp1

        den = 1.0 - x0 * x0
        d1 = n * (pkm1 - x0 * pk)
        dpn = d1 / den
        d2pn = (2.0 * x0 * dpn - e1 * pk) / den
        d3pn = (4.0 * x0 * d2pn + (2.0 - e1) * dpn) / den
        d4pn = (6.0 * x0 * d3pn + (6.0 - e1) * d2pn) / den
        u = pk / dpn
        v = d2pn / dpn
        h = -u * (1.0 + 0.5 * u * (v + u * (v * v - d3pn / (3.0 * dpn))))
        p = pk + h * (dpn + 0.5 * h * (d2pn + h / 3.0 * (d3pn + 0.25 * h * d4pn)))
        dp = dpn + h * (d2pn + 0.5 * h * (d3pn + h * d4pn / 3.0))
        h = h - p / dp
        fx = d1 - h * e1 * (pk + 0.5 * h * (dpn + h / 3.0 * (d2pn + 0.25 * h * (d3pn + 0.2 * h * d4pn))))
        x.append(x0 + h)
        w.append(2.0 * (1.0 - x[i] * x[i]) / (fx * fx))

    if m + m > n:
        x[m - 1] = 0.0

    return x, w

def main_solution(n):
    x, w = grule(n)
    x_serializable = [float(val) for val in x]
    w_serializable = [float(val) for val in w]
    return {"roots": x_serializable, "weights": w_serializable}

# Test the function with n = 9
result = main_solution(9)
print(result)