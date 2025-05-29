from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 92 - 94*x
poly2 = 31*x + 59
poly3 = -32*x**3 - 64*x + 72

# Expand the first two polynomials
expanded_first_two = expand(poly1 * poly2)

# Expand the result with the third polynomial
final_expansion = expand(expanded_first_two * poly3)

# Print the final expanded polynomial
print(final_expansion)