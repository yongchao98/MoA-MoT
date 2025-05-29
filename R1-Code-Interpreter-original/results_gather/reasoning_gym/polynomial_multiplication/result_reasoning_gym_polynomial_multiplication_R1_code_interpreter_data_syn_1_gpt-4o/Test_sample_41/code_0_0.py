from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = 54*y**3 - 64*y
poly2 = -3*y**3 - 129*y**2 - 76

# Perform the multiplication and simplify
result = expand(poly1 * poly2)

# Output the result
print(result)