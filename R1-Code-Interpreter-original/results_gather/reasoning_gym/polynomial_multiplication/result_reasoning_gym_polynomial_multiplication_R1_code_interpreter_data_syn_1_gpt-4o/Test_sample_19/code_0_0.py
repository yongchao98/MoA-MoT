from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 51*x**3 + 30*x**2
poly2 = 47*x**3 + 64*x**2 + 52

# Perform the multiplication and simplify
result = expand(poly1 * poly2)

# Print the simplified result
print(result)