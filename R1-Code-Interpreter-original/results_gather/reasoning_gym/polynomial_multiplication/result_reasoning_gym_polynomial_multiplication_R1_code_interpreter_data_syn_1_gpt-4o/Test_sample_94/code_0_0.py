from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = -38*z**2 + 18*z
poly2 = 30*z**2 - 135*z

# Perform the multiplication and simplify
result = expand(poly1 * poly2)

# Print the result
print(result)