from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 25 - 87*x**2
poly2 = 14*x**2 + 50*x + 89

# Expand the product of the two polynomials
expanded_poly = expand(poly1 * poly2)

# Print the expanded polynomial
print(expanded_poly)