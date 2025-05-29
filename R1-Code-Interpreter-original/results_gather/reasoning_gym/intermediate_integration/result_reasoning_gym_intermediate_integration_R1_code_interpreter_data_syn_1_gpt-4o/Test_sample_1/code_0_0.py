from sympy import symbols, integrate, sin, cos

x = symbols('x')
integrand = 3*x**2*sin(3*x)
integral = integrate(integrand, x)
print(integral)