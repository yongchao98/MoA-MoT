from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expression
expression = (29 - 49*y**2)*(21*y**3 + 49*y)

# Expand the expression
expanded_expression = expand(expression)

# Print the expanded expression
print(expanded_expression)