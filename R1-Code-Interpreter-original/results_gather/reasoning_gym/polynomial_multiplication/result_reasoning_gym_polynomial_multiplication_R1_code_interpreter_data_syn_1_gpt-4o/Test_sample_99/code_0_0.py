from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = 17*y + 27
poly2 = 213*y**2 - 95
poly3 = -77*y**3 + 34*y - 20

# Expand the expression
expanded_expr = expand(poly1 * poly2 * poly3)

# Print the expanded expression
print(expanded_expr)