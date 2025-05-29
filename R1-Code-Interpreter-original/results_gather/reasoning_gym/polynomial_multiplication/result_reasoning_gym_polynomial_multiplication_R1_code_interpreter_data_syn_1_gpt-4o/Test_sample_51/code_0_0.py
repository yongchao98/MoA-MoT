from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = -84*y**2 - 17*y
poly2 = 74*y**2 - 67*y - 6
poly3 = 46*y**3 + 26*y - 11

# Multiply the first two polynomials
intermediate_poly = expand(poly1 * poly2)

# Multiply the result with the third polynomial
final_result = expand(intermediate_poly * poly3)

# Print the final expanded polynomial
print(final_result)