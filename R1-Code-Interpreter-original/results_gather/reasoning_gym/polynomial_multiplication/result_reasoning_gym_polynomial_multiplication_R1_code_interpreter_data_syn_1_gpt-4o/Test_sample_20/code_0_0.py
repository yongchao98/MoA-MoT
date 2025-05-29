from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expression
expr = (71*y**2 + 18*y)*(-9*y**2 - 50*y - 81)

# Expand the expression
expanded_expr = expand(expr)

# Print the expanded expression
print(expanded_expr)