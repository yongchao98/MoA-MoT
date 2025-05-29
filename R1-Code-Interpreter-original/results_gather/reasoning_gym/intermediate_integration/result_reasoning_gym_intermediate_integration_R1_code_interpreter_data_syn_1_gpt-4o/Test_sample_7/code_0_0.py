from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = 12 * sqrt(4*x + 7)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)