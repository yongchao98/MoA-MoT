from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expressions
expr1 = -81*y**2 - 37
expr2 = -51*y**2 - 147*y + 35

# Expand the product of the two expressions
simplified_expr = expand(expr1 * expr2)

# Print the simplified expression
print(simplified_expr)