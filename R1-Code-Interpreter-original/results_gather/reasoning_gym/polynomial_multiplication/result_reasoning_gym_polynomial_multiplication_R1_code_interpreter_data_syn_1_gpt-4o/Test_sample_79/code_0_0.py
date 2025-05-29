from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = -58*x**3 + 41*x
poly2 = -25*x**3 - 50*x
poly3 = 54*x**3 - 93*x**2 + 49*x

# First multiplication
result1 = expand(poly1 * poly2)

# Second multiplication
final_result = expand(result1 * poly3)

# Print the final expanded polynomial
print(final_result)