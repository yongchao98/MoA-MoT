from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the expression
expr = (76*x - 20)*(18*x**2 - 11)*(44*x**3 - 62*x**2 - 150)

# Expand the expression
expanded_expr = expand(expr)

# Print the expanded expression
print(expanded_expr)