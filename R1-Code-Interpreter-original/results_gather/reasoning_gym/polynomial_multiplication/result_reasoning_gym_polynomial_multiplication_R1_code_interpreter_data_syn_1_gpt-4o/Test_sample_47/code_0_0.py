from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = -112*x**3 - x**2
poly2 = -93*x**3 - 77*x - 45
poly3 = 35*x**3 - 35*x**2 - 15*x

# Multiply the polynomials
result = expand(poly1 * poly2 * poly3)

# Print the result
print(result)