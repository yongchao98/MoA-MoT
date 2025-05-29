from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 107*x**3 + 16
poly2 = 88*x**3 + 99*x**2 - 25*x

# Perform the multiplication
result = expand(poly1 * poly2)

# Print the result
print(result)