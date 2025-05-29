from sympy import symbols, integrate, cos, sin

# Define the variable
x = symbols('x')

# Define the function to integrate
f = -x * cos(2*x)

# Perform the integration
integral_result = integrate(f, x)

# Print the result
print(integral_result)