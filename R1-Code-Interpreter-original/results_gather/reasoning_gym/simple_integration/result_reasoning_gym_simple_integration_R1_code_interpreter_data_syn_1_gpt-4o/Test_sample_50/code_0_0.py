from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the expression
expression = 45*X**8 + 56*X**7 + 15*X**4/4

# Integrate the expression
integral = integrate(expression, X)

# Print the result
print(integral)