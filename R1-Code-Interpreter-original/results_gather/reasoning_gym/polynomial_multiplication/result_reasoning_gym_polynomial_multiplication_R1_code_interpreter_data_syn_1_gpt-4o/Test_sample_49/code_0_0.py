from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = -65*z**3 + 88*z**2
poly2 = -3*z**2 - 72*z + 87
poly3 = -24*z**3 + 157*z**2 - 3*z

# Multiply the first two polynomials
result1 = expand(poly1 * poly2)

# Multiply the result with the third polynomial
final_result = expand(result1 * poly3)

# Print the final expanded polynomial
print(final_result)