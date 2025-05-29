from sympy import symbols, expand

x = symbols('x')
expression = (-107*x**3 - 120*x)*(12*x**3 + 50)
expanded_expression = expand(expression)
print(expanded_expression)