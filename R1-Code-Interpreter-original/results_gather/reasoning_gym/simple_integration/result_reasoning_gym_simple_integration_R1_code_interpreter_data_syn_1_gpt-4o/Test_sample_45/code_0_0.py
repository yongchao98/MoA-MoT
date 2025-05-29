from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function to integrate
function = 10*X**9 + 81*X**8 + 16*X**3 + 12*X**2

# Perform the integration
integral_result = integrate(function, X)

# Print the result
print(integral_result)