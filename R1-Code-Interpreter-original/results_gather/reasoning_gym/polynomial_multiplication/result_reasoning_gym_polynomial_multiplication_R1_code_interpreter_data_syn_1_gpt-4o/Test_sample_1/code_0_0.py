from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = 36*z**3 + 136*z - 82
poly2 = 49*z**3 + 69*z**2 - 68*z

# Perform the multiplication and simplify
result = expand(poly1 * poly2)

# Print the result
print(result)