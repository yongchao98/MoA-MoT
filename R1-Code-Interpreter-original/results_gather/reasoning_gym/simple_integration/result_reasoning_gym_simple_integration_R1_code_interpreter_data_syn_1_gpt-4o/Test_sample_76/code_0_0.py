from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 154*x**3/15

# Perform the integration
integral_result = integrate(function, x)

# Print the result
print(integral_result)