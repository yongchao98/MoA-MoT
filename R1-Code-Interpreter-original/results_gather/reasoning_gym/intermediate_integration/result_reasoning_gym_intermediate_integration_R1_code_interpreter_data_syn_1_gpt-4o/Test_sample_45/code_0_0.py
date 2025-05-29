from sympy import symbols, exp, sin, cos, integrate

x = symbols('x')
integral = integrate(exp(x)*sin(x), x)
print(integral)