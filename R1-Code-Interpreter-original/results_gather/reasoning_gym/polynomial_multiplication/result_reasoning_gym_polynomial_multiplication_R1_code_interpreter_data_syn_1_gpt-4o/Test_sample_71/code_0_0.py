from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the polynomials
poly1 = -85*z**2 - 51*z
poly2 = 14*z**2 - 75
poly3 = 60*z**3 + 90*z**2 - 16*z + 78

# Perform the multiplications
result = expand(poly1 * poly2 * poly3)

# Print the simplified result
print(result)