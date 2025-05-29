from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the expression
expression = 80*X**9 + 9*X**8 - 376*X**7/5 + 21*X**6/5

# Integrate the expression
integral = integrate(expression, X)

# Print the result
print(integral)