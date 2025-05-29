from sympy import symbols, exp, integrate

x = symbols('x')
integral = integrate(2*x**3*exp(x), x)
print(integral)