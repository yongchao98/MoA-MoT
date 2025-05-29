from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = 71*z**2 + 29*z + 51
poly2 = -3*z**3 + 44*z**2 + 95*z
poly3 = -97*z**3 + 84*z**2 - 3*z + 68

# Multiply the first two polynomials
result1 = expand(poly1 * poly2)

# Multiply the result with the third polynomial
final_result = expand(result1 * poly3)

# Print the simplified expression
print(final_result)