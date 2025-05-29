from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function to integrate
function = 18*X**8 - 15*X**2

# Perform the integration
integral_result = integrate(function, X)

# Print the result
print(integral_result)