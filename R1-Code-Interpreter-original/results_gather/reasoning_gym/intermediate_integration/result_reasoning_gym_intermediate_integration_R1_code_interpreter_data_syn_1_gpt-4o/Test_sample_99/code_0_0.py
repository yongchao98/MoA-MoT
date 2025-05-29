from sympy import symbols, integrate, sin, cos

# Define the variable and constant
x = symbols('x')
a = symbols('a')

# Define the function to integrate
function = -a * sin(2*x)

# Perform the integration
integral_result = integrate(function, x)

# Print the result
print(integral_result)