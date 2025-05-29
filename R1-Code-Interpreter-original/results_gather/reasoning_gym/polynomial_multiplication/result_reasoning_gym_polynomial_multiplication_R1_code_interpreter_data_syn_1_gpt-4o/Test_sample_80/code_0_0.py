from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = 50 - 87*y**2
poly2 = -67*y**3 + 112*y**2
poly3 = -21*y**3 + 29*y - 59

# Multiply the first two polynomials
result1 = expand(poly1 * poly2)

# Multiply the result with the third polynomial
final_result = expand(result1 * poly3)

# Print the simplified expression
print(final_result)