from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = -16*z**3 + 4*z
poly2 = 49*z**3 + 117*z**2 + 88

# Expand the product of the two polynomials
simplified_expression = expand(poly1 * poly2)

# Print the simplified expression
print(simplified_expression)