from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = -72*z - 34
poly2 = 84*z**3 - 143*z
poly3 = -37*z**3 - 83*z**2 - 47*z

# Multiply the first two polynomials
result1 = expand(poly1 * poly2)

# Multiply the result with the third polynomial
final_result = expand(result1 * poly3)

# Print the simplified expression
print(final_result)