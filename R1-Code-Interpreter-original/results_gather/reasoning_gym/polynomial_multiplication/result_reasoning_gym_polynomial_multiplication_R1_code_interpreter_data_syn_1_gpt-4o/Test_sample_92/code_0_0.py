from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expression
expression = (-101*y - 28)*(-8*y - 47)

# Simplify the expression
simplified_expression = expand(expression)

# Print the simplified expression
print(simplified_expression)