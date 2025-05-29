from sympy import symbols, expand

x = symbols('x')
expression = 53*x**3*(-16*x**3 + 83*x)*(-98*x**3 + 95*x - 78)
simplified_expression = expand(expression)
print(simplified_expression)