from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = -43*z**2 + 89*z
poly2 = 115*z**2 + 18*z
poly3 = 10*z**3 + 75*z**2 - 12*z

# Expand the product of the first two polynomials
expanded_poly1_2 = expand(poly1 * poly2)

# Expand the result with the third polynomial
expanded_expression = expand(expanded_poly1_2 * poly3)

# Print the simplified expression
print(expanded_expression)