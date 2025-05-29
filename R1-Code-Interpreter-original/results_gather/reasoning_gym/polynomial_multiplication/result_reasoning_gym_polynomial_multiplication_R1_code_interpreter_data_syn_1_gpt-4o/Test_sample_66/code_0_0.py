from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expression
expression = (91 - 26*y) * (-197*y**3 - 83*y + 100)

# Expand the expression
simplified_expression = expand(expression)

# Print the simplified expression
print(simplified_expression)