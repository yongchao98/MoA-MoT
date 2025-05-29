from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = -14*x**3 - 38
poly2 = -9*x**3 - 51*x**2
poly3 = 33*x**2 - 4*x + 94

# Expand the product of the first two polynomials
result1 = expand(poly1 * poly2)

# Expand the product of the result with the third polynomial
final_result = expand(result1 * poly3)

# Print the final expanded polynomial
print(final_result)