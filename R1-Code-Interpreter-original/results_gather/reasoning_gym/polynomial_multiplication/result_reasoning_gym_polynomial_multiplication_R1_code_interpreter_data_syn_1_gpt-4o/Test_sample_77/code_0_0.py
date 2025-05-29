from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expression
expr = (65*y + 15)*(-105*y**3 + 86*y**2)

# Expand the expression
simplified_expr = expand(expr)

# Print the simplified expression
print(simplified_expr)