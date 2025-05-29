from sympy import symbols, integrate, sin, cos

x = symbols('x')
integral = integrate(x**3 * sin(2*x), x)
print(integral)