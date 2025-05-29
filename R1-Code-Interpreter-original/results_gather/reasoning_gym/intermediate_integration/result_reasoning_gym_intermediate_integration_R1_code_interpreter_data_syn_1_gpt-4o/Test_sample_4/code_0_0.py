from sympy import symbols, integrate, asin

x = symbols('x')
integral = integrate(asin(x), x)
print(integral)