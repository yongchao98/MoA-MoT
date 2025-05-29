from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the expression
expression = 99 * x**3 * (8 * x**3 - 68 * x**2 - 84)

# Expand the expression
expanded_expression = expand(expression)

# Print the expanded expression
print(expanded_expression)