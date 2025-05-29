from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 53*x**2 - 46
poly2 = -26*x**3 + 52*x**2
poly3 = -13*x**2 - 47*x - 15

# Perform the multiplications
result = expand(poly1 * poly2 * poly3)

# Print the expanded and simplified result
print(result)