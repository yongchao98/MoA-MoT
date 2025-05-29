from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = -137*x**3 + 36*x
poly2 = -13*x**3 - 23*x**2 - 57*x - 52

# Perform the multiplication and simplify
result = expand(poly1 * poly2)

# Print the result
print(result)