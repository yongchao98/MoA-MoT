from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 86*x**2 + 57*x
poly2 = 44*x**3 - 31*x + 49

# Perform the multiplication and simplify
result = expand(poly1 * poly2)

# Print the result
print(result)