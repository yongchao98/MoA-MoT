from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = -34*y - 9
poly2 = -15*y**3 + 31*y**2
poly3 = -79*y**2 - 29*y + 81

# Perform the multiplication
result = expand(poly1 * poly2 * poly3)

# Print the simplified expression
print(result)