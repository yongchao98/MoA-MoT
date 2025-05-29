from sympy import symbols, sin, cos, integrate

# Define the variable
x = symbols('x')

# Define the function to integrate
integrand = -12 * sin(6*x - 8) * cos(6*x - 8)**3

# Perform the integration with respect to x
antiderivative = integrate(integrand, x)

# Display the result
print(antiderivative)