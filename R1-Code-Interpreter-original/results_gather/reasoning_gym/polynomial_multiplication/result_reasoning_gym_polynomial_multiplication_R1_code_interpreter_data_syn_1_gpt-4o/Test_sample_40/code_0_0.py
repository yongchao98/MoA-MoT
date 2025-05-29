from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the expression
expression = (-180*x**2 + 92*x)*(29*x**3 + 82)

# Expand the expression
expanded_expression = expand(expression)

# Print the expanded expression
print(expanded_expression)