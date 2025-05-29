from sympy import symbols, integrate, asin, sqrt

# Define the variable and function
x = symbols('x')
f = -2 * asin(x)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)