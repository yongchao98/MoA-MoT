from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the expression
expression = 90*X**9 + 56*X**7 + 35*X**6/3

# Calculate the indefinite integral
indefinite_integral = integrate(expression, X)

# Print the result
print(indefinite_integral)