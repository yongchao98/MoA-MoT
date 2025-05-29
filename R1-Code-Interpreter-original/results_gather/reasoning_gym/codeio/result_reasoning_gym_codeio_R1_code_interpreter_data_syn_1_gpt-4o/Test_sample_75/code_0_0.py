import math

def my_func(x):
    return 10 * math.cos(x) - 0.1 * x * x

def roots_step(a, b, h, func):
    arr = []
    l = a
    r = l + h
    while l < b:
        y_1 = func(l)
        y_2 = func(r)
        if y_1 * y_2 <= 0:
            arr.append((l, r))
        l = r
        r = r + h
        if r > b:
            r = b
    return arr

def method_bisection(interval_root, func, epsilon=1e-5):
    root = []
    for i in interval_root:
        a = i[0]
        b = i[1]
        while (b - a) > epsilon:
            c = (a + b) / 2
            if func(a) * func(c) <= 0:
                b = c
            else:
                a = c
        X = (a + b) / 2
        root.append(X)
    return root

def main_solution(a, b, h, epsilon):
    func = my_func
    intervals = roots_step(a, b, h, func)
    roots = method_bisection(intervals, func, epsilon)
    return roots

# Test different intervals and step sizes
a = 4.0
b = 8.0
h = 0.5
epsilon = 1e-5

roots = main_solution(a, b, h, epsilon)
print(roots)