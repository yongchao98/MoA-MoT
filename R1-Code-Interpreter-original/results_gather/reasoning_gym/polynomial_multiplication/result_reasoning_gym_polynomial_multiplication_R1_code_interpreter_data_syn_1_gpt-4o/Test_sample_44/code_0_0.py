from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
expr1 = (32 - 155*y)
expr2 = (-57*y**3 - 57*y)
expr3 = (-39*y**2 + 53*y - 69)

# Expand the expression
expanded_expr = expand(expr1 * expr2 * expr3)

# Print the expanded expression
print(expanded_expr)