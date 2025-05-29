from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the expression
expression = 90*X**8 + 12*X**5 - 24*X**2 + 29*X/2

# Calculate the indefinite integral
indefinite_integral = integrate(expression, X)

# Print the result
print(indefinite_integral)