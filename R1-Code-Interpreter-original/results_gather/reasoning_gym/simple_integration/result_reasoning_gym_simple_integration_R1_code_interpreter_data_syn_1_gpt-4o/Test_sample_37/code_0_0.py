from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function to integrate
function = -20*x**9 - 28*x**7/5

# Perform the integration
integral_result = integrate(function, x)

# Print the result
print(integral_result)