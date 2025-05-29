from sympy import symbols, integrate, cos, sin

x = symbols('x')
integral = integrate(-2*x**3*cos(x), x)
print(integral)