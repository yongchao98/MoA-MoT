from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 4*x**2 - 53*x
poly2 = 40*x**3 + 41*x**2 + 140

# Perform the multiplication
result = expand(poly1 * poly2)

# Print the result
print(result)