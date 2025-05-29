from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = 71 - 8*z**3
poly2 = -13*z**3 + z**2
poly3 = 50*z**3 + 4*z**2 + 9*z + 44

# Perform the multiplication and expansion
result = expand(poly1 * poly2 * poly3)

# Print the simplified expression
print(result)