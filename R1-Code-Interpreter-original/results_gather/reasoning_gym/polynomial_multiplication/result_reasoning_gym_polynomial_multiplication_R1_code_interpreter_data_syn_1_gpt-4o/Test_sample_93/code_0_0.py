from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 24*x**3 - 38*x + 82
poly2 = 20*x**3 + 78*x**2 - 39*x + 1

# Expand the product of the two polynomials
expanded_poly = expand(poly1 * poly2)

# Print the expanded polynomial
print(expanded_poly)