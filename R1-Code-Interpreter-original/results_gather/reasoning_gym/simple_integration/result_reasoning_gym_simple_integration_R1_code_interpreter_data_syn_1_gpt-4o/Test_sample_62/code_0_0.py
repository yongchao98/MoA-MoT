from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the expression
expression = -35*x**9/3 - 49/6

# Calculate the antiderivative
antiderivative = integrate(expression, x)

# Print the result
print(antiderivative)