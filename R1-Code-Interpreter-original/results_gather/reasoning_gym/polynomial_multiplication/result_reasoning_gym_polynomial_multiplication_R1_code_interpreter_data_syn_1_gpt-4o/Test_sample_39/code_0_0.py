from sympy import symbols, expand

y = symbols('y')
expression = (91 - 20*y**2)*(30*y**2 + 47*y - 78)
expanded_expression = expand(expression)
print(expanded_expression)