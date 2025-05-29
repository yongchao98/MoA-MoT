from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = 35*z**3 - 58*z
poly2 = 78*z**2 - 90*z - 97
poly3 = -14*z**3 - 68*z**2 - 11*z - 48

# Expand the expression
expanded_expr = expand(poly1 * poly2 * poly3)

# Print the simplified expression
print(expanded_expr)