from sympy import symbols, expand

y = symbols('y')
expression = 12*y*(-66*y**2 + 112*y)*(95*y**2 - 33)
simplified_expression = expand(expression)
print(simplified_expression)