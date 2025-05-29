from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = 49*z**2 - z
poly2 = -78*z**2 - 94*z - 144
poly3 = -164*z**3 - 32*z + 96

# Multiply the first two polynomials
result1 = expand(poly1 * poly2)

# Multiply the result with the third polynomial
final_result = expand(result1 * poly3)

# Print the simplified expression
print(final_result)