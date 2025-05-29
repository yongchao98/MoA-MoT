from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = 72*z**3 - 158
poly2 = 96*z**3 + 37*z
poly3 = -31*z**3 + 44*z**2 - 87*z

# Multiply the first two polynomials
result1 = expand(poly1 * poly2)

# Multiply the result with the third polynomial
final_result = expand(result1 * poly3)

# Print the final expanded polynomial
print(final_result)