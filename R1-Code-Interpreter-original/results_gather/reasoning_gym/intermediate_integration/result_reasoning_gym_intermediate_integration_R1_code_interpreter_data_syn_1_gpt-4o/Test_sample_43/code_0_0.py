from sympy import symbols, integrate, sin, cos

# Define the variable and the function
x = symbols('x')
f = 27 * sin(9*x - 6) * cos(9*x - 6)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)