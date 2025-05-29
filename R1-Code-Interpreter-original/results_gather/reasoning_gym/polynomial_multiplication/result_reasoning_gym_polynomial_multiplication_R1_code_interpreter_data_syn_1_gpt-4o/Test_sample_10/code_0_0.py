from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the expression
expr = (18 - 95*z**3)*(-12*z**2 - 104*z + 78)

# Expand the expression
expanded_expr = expand(expr)

# Print the expanded expression
print(expanded_expr)