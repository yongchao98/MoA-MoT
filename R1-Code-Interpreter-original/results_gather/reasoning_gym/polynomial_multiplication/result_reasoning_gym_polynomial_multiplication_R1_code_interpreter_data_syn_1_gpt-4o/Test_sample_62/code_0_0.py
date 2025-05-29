from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the expression
expr = (-36*z**3 - 3*z)*(111*z**3 - 3*z)

# Expand the expression
expanded_expr = expand(expr)

# Print the expanded expression
print(expanded_expr)