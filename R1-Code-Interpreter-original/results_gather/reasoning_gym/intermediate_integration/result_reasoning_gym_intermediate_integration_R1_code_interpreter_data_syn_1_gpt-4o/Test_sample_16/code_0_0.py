from sympy import symbols, integrate, cos, sin

# Define the variable
x = symbols('x')

# Define the integrand
integrand = 2*x**3*cos(3*x)

# Perform the integration
integral_result = integrate(integrand, x)

# Print the result
print(integral_result)