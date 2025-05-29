from sympy import symbols, cos, integrate

# Define the variable and function
x = symbols('x')
u = cos(6*x - 8)

# Perform the integration
antiderivative = integrate(2 * u**3, u)

# Display the result
print(antiderivative)