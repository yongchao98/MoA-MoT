from sympy import symbols, expand

z = symbols('z')
expression = (39 - 126*z**2)*(93*z**2 - 20*z)
expanded_expression = expand(expression)
print(expanded_expression)