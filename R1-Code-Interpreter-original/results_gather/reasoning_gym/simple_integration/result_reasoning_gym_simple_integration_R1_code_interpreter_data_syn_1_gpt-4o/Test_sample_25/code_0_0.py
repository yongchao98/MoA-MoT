from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = -36*X**8 + 49*X**6/6

# Perform the integration
integral_result = integrate(function, X)

# Print the result
print(integral_result)