from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 74*x**2 + 98*x + 16
poly2 = 77*x**2 + 8*x - 49
poly3 = -163*x**3 - 49*x + 8

# Expand the expression
expanded_expr = expand(poly1 * poly2 * poly3)

# Print the simplified expression
print(expanded_expr)