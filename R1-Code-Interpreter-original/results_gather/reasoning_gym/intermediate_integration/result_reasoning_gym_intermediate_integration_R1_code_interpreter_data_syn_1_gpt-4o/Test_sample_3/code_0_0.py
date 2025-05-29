from sympy import symbols, integrate, sin, cos

x = symbols('x')
integrand = -2*x**3*sin(x)
antiderivative = integrate(integrand, x)
print(antiderivative)