from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
poly1 = 31*x**2 + 63*x
poly2 = -16*x**3 + 42*x**2
poly3 = -5*x**3 - 88*x**2 + 25*x

# Multiply the first two polynomials
result1 = expand(poly1 * poly2)

# Multiply the result by the third polynomial
final_result = expand(result1 * poly3)

# Print the simplified expression
print(final_result)