from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = -49*x**3 + 93*x
poly2 = 78*x**3 + 11
poly3 = 80*x**3 - 71*x

# Expand the product of the three polynomials
expanded_expression = expand(poly1 * poly2 * poly3)

# Print the expanded expression
print(expanded_expression)