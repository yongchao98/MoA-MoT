from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = 10 * sqrt(10*x + 9)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)