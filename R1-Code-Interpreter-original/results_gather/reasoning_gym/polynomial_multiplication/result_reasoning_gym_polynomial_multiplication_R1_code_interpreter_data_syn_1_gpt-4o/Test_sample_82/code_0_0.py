from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 68*x**2 + 93*x
poly2 = 56*x**3 + 125*x**2

# Perform the multiplication and simplify
result = expand(poly1 * poly2)

# Print the result
print(result)