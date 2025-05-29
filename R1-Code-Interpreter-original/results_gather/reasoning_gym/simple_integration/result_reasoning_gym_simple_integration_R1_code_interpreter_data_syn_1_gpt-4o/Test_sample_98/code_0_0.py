from sympy import symbols, integrate

X = symbols('X')
f = (54/5)*X**8 - 6*X**7 + (2/3)*X**5 + 18*X**2 - 6
antiderivative = integrate(f, X)
print(antiderivative)