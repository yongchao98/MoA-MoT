from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = 62*z + 97
poly2 = 10*z**2 + 11*z + 67
poly3 = -41*z**3 + 53*z - 74

# Expand the first two polynomials
expanded_poly1_2 = expand(poly1 * poly2)

# Expand the result with the third polynomial
final_expansion = expand(expanded_poly1_2 * poly3)

# Print the final expanded polynomial
print(final_expansion)