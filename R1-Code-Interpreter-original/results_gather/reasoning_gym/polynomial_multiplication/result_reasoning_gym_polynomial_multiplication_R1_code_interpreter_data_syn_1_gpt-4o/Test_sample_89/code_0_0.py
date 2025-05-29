from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the expression
expr = (155*x + 55)*(-55*x**2 - 86*x - 56)

# Expand the expression
simplified_expr = expand(expr)

# Print the simplified expression
print(simplified_expr)