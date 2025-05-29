from sympy import symbols, log, diff

# Define the variable and function
x = symbols('x')
integral_result = 3*x*log(x**2) - 6*x

# Differentiate the result
derivative = diff(integral_result, x)

# Output the derivative
print(derivative)