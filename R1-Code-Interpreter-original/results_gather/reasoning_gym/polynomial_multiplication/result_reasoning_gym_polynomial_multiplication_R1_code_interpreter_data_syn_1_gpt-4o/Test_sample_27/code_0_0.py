from sympy import symbols, expand

x = symbols('x')
expression = (75*x**2 + 20)*(-2*x**3 + 69*x**2 - 20)
simplified_expression = expand(expression)
print(simplified_expression)