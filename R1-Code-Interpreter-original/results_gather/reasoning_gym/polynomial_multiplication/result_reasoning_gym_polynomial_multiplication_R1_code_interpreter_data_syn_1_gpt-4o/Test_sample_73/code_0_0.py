from sympy import symbols, expand

# Define the variable
z = symbols('z')

# Define the expression
expression = 11 * z**3 * (87 * z**2 + 81 * z) * (-99 * z**3 + 45 * z**2 + 5)

# Expand the expression
expanded_expression = expand(expression)

# Print the expanded expression
print(expanded_expression)