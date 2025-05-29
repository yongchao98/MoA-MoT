from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = -65*y**2 + 79*y - 88
poly2 = 57*y**2 + 11*y - 1

# Expand the product of the two polynomials
expanded_expression = expand(poly1 * poly2)

# Print the expanded expression
print(expanded_expression)