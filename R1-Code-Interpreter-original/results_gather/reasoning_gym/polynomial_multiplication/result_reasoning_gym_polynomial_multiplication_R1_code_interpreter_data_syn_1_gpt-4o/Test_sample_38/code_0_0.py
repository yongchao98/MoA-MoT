from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 69*x - 188
poly2 = -112*x**3 + 101*x
poly3 = 49*x**3 + 10*x**2 - 37*x - 12

# Perform the multiplication and expansion
result = expand(poly1 * poly2 * poly3)

# Print the simplified result
print(result)