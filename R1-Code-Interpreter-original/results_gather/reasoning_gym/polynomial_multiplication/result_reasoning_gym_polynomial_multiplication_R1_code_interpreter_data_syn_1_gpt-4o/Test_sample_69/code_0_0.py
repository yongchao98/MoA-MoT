from sympy import symbols, expand

z = symbols('z')
expr = (-57*z**3 - 33*z**2)*(28*z**3 + 65*z**2 - 9*z)
expanded_expr = expand(expr)
print(expanded_expr)