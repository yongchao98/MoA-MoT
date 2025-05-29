from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = 51 - 12*y**3
poly2 = 84*y**2 + 48*y + 89
poly3 = -77*y**3 + 62*y**2 + 99*y - 9

# Perform the multiplication and expansion
result = expand(poly1 * poly2 * poly3)

# Print the simplified result
print(result)