from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 59 - 106*x
poly2 = 95*x**2 + 116
poly3 = -61*x**3 + 62*x

# Expand the first two polynomials
expanded_poly1_2 = expand(poly1 * poly2)

# Multiply the result with the third polynomial
final_expansion = expand(expanded_poly1_2 * poly3)

# Print the final expanded polynomial
print(final_expansion)