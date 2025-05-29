from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expressions
expr1 = 20*y**2 + 9
expr2 = 22*y**3 - 156
expr3 = 89*y**3 + 22*y

# Multiply the first two expressions
result1 = expand(expr1 * expr2)

# Multiply the result with the third expression
final_result = expand(result1 * expr3)

# Print the final expanded expression
print(final_result)