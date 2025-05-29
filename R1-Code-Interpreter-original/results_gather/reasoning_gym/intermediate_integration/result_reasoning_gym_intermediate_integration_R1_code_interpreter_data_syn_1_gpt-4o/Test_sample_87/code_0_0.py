from sympy import symbols, integrate, sin, cos

x = symbols('x')
integral = integrate(-x**3 * sin(x), x)
print(integral)