from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expression
expression = (3*y**2 + 26*y) * (45*y**3 - 37*y**2)

# Expand the expression
expanded_expression = expand(expression)

# Print the expanded expression
print(expanded_expression)