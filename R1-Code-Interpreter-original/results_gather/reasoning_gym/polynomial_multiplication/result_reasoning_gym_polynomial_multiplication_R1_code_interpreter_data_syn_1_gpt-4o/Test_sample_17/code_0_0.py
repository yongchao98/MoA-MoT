from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 66 - 83*x**2
poly2 = -42*x**2 + 89*x + 52
poly3 = 78*x**2 + 6*x - 70

# Expand the product of the three polynomials
expanded_expr = expand(poly1 * poly2 * poly3)

# Print the simplified expression
print(expanded_expr)