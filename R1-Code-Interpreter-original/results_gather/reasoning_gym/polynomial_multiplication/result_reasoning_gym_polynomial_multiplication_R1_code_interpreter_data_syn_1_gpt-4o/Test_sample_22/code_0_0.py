from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = -2*x**2 - 106
poly2 = 19*x**3 - 27*x
poly3 = 150*x**2 - 35*x + 58

# Expand the product of the first two polynomials
expanded_poly1_2 = expand(poly1 * poly2)

# Expand the result with the third polynomial
expanded_full = expand(expanded_poly1_2 * poly3)

# Print the simplified expression
print(expanded_full)