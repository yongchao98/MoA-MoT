from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the expression
expr = (95*x**2 - 96*x)*(43*x**3 + 55)

# Expand the expression
simplified_expr = expand(expr)

# Print the simplified expression
print(simplified_expr)