from sympy import symbols, exp, integrate

x = symbols('x')
integrand = (12 - 54*x) * exp(9*x**2 - 4*x - 8)
integral = integrate(integrand, x)
print(integral)