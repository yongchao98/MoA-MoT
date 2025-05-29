from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = 21*sqrt(7*x - 7)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)