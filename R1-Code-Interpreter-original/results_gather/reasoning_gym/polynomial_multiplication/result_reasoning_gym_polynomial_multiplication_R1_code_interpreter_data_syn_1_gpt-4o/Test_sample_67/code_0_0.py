from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = -77*y**2 - 18
poly2 = 16*y**2 + 68
poly3 = 35*y**2 + 44

# Expand the product of the first two polynomials
expanded_poly1_2 = expand(poly1 * poly2)

# Expand the result with the third polynomial
final_expansion = expand(expanded_poly1_2 * poly3)

# Print the final expanded expression
print(final_expansion)