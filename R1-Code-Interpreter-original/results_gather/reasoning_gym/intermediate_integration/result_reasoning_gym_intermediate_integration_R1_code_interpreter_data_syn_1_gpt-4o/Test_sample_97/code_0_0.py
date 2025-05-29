from sympy import symbols, integrate, cos, sin

x = symbols('x')
integral = integrate(3*x*cos(2*x), x)
print(integral)