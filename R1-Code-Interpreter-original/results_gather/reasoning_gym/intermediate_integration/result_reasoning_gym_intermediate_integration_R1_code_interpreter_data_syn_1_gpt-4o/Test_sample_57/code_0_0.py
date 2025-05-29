from sympy import symbols, integrate

# Define the variable and the expression
x = symbols('x')
expression = -2 * (3*x + 6)**2

# Perform the integration
integral_result = integrate(expression, x)

# Print the result
print(integral_result)