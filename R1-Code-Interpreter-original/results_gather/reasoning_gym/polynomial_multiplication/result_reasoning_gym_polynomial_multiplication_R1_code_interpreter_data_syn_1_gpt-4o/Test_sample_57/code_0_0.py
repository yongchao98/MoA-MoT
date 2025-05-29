from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 54*x**2 - 45
poly2 = 72*x**2 - 12
poly3 = 22*x**3 - 89*x**2 + 29*x

# Multiply the first two polynomials
result1 = expand(poly1 * poly2)

# Multiply the result by the third polynomial
final_result = expand(result1 * poly3)

# Print the final simplified expression
print(final_result)