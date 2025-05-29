from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expression
expression = -25*y**3 * (-99*y**3 + 123*y**2) * (9*y**3 + 43*y)

# Expand the expression
expanded_expression = expand(expression)

# Print the expanded expression
print(expanded_expression)