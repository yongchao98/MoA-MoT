from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the polynomials
expr1 = 74 - 123*y**3
expr2 = 96 - 104*y**2
expr3 = 15*y**3 - 36

# Expand the first two polynomials
expanded_expr1_2 = expand(expr1 * expr2)

# Expand the result with the third polynomial
final_expansion = expand(expanded_expr1_2 * expr3)

# Print the simplified expression
print(final_expansion)