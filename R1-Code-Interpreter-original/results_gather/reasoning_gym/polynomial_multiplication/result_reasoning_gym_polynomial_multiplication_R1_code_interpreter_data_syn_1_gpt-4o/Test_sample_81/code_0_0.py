from sympy import symbols, expand

# Define the variable
x = symbols('x')

# Define the polynomials
expr1 = (40 - 25*x**3)
expr2 = (-62*x**2 - 4)
expr3 = (83*x**3 + 80*x)

# Perform the multiplications
result = expand(expr1 * expr2 * expr3)

# Print the simplified expression
print(result)