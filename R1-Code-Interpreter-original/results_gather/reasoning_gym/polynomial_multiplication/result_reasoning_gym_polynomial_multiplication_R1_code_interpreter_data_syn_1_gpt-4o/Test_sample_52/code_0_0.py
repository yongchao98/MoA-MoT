from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the expression
expression = 186 * z**3 * (-6 * z**2 - 85 * z) * (85 * z**3 + 63)

# Expand the expression
expanded_expression = expand(expression)

# Print the expanded expression
print(expanded_expression)