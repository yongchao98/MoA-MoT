from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = -93*y**2 - 82
poly2 = -14*y**3 + 8*y**2 + 72
poly3 = -20*y**3 - 24*y**2 + 62*y + 6

# Expand the expression
expanded_expr = expand(poly1 * poly2 * poly3)

# Print the simplified expression
print(expanded_expr)