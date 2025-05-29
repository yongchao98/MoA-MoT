from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the expression
expression = 70*X**6/3 + 10*X**3 + 24*X**2

# Compute the indefinite integral
indefinite_integral = integrate(expression, X)

# Print the result
print(indefinite_integral)