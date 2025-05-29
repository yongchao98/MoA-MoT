from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the expression
expression = -103*z**3*(6*z**2 - 30*z + 36)

# Expand the expression
expanded_expression = expand(expression)

# Print the expanded expression
print(expanded_expression)