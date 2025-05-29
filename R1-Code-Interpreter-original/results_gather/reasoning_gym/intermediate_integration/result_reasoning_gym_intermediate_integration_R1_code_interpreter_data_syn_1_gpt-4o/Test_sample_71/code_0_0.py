from sympy import symbols, integrate, exp

x = symbols('x')
integral = integrate(x**3 * exp(x), x)
print(integral)