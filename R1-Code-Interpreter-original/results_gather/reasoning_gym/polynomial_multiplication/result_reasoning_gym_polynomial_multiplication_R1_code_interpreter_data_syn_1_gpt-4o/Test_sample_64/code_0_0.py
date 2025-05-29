from sympy import symbols, expand

# Define the variable
y = symbols('y')

# Define the expressions
expr1 = (3 - 54*y**2)
expr2 = (10 - 73*y**3)
expr3 = (-46*y**3 - 60*y)

# Multiply and expand the expressions
result = expand(expr1 * expr2 * expr3)

# Print the result
print(result)