from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -(2*x + 6)**4

# Perform the integration
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)