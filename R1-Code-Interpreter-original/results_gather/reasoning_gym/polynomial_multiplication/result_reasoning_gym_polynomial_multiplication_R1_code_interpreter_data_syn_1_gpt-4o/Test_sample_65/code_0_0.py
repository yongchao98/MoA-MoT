from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the expression
expr = -127 * x * (13 - 35 * x**2) * (57 * x**3 + 24 * x**2 - 28 * x)

# Expand the expression
expanded_expr = expand(expr)

# Print the expanded expression
print(expanded_expr)