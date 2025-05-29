from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
expr1 = -89*y - 80
expr2 = 80*y**2 - 68
expr3 = 48*y**3 - 81*y**2

# Expand the expression
expanded_expr = expand(expr1 * expr2 * expr3)

# Print the simplified expression
print(expanded_expr)