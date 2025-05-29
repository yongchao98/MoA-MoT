from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 34 - 43*x
poly2 = -25*x**2 - 53*x
poly3 = -94*x**2 + 99*x - 15

# First multiplication
result1 = expand(poly1 * poly2)

# Second multiplication
final_result = expand(result1 * poly3)

# Print the final expanded polynomial
print(final_result)