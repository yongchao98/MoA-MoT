from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function to integrate
function = 2 * (x - 8)**2

# Perform the integration
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)