from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function to integrate
function = -8*X**4 - 24*X**2

# Perform the integration
integral_result = integrate(function, X)

# Print the result
print(integral_result)