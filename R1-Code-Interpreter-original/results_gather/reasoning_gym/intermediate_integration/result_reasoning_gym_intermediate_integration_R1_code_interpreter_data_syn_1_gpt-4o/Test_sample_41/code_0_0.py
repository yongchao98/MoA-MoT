from sympy import symbols, atan, integrate, ln

# Define the variable
x = symbols('x')

# Define the function to integrate
f = -2 * atan(x)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)