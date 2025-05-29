from sympy import symbols, expand

x = symbols('x')
expression = (-64*x**2 - 11)*(81*x**2 + 50)
expanded_expression = expand(expression)
print(expanded_expression)