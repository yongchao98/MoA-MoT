from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expression
expression = 140*y**2 * (y**2 - 67*y) * (-10*y**3 + 65*y**2 + 98)

# Expand the expression
expanded_expression = expand(expression)

# Print the expanded expression
print(expanded_expression)