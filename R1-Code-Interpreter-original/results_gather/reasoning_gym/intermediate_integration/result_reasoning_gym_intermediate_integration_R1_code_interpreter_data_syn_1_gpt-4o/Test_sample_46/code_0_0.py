from sympy import symbols, exp, cos, sin, integrate

x = symbols('x')
integral = integrate(-3*exp(x)*cos(x), x)
print(integral)