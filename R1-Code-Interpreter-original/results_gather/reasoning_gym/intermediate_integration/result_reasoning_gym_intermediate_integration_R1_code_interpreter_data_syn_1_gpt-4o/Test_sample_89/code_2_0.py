from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the expression
u = 3*x + 10

# Expand the expression (3*x + 10)**4
expanded_expression = expand(u**4)

# Print the expanded expression
print(expanded_expression)