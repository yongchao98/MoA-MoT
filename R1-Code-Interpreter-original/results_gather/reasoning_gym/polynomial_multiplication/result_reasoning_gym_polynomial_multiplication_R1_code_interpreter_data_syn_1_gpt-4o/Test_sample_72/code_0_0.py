from sympy import symbols, expand

x = symbols('x')
expression = (-101*x - 19) * (-75*x**2 + 42*x + 48)
expanded_expression = expand(expression)
print(expanded_expression)