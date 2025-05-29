from sympy import symbols, integrate, cos, sin

x = symbols('x')
integral = integrate(-2*x**2*cos(3*x), x)
print(integral)