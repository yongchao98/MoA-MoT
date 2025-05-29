from sympy import symbols, sin, cos, integrate

# Define the variable and constant
x = symbols('x')
a = symbols('a')

# Define the function to integrate
function = -a * sin(2 * x)

# Perform the integration
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)