from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = -26*y**3 + 33*y - 50
poly2 = -17*y**3 - 77*y**2 + 37
poly3 = 44*y**3 - 25*y - 14

# Expand the product of the polynomials
expanded_poly = expand(poly1 * poly2 * poly3)

# Print the expanded polynomial
print(expanded_poly)