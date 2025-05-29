from sympy import symbols, integrate, cos, sin

x = symbols('x')
integral = integrate(-x**3 * cos(x), x)
print(integral)