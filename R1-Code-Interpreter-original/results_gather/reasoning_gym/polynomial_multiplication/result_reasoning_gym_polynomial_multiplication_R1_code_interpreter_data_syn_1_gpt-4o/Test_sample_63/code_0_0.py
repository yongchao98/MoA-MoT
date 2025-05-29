from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = 4*z**2 + 102
poly2 = -64*z**3 + 78*z**2 + 24*z

# Perform the multiplication and simplify
result = expand(poly1 * poly2)

# Print the simplified result
print(result)