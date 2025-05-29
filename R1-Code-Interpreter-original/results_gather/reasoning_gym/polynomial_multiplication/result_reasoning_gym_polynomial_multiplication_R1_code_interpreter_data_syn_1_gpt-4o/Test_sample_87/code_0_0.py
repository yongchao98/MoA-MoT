from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
poly1 = -44*y**2 + 2*y + 56
poly2 = -33*y**3 + 55*y**2 - 69*y - 75

# Perform the multiplication and simplify
result = expand(poly1 * poly2)

# Print the result
print(result)