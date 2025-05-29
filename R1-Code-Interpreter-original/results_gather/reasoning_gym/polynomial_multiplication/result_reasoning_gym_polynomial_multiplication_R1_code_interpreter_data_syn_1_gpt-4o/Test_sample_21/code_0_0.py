from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = 22*z**3 - 34
poly2 = 67*z**3 + 65*z + 94

# Perform the multiplication
result = expand(poly1 * poly2)

# Print the result
print(result)