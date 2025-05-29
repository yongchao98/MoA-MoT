from sympy import symbols, integrate, sqrt

# Define the variable and the function
x = symbols('x')
f = -8 * sqrt(8 - 4*x)

# Perform the integration
integral_result = integrate(f, x)

# Output the result
print(integral_result)