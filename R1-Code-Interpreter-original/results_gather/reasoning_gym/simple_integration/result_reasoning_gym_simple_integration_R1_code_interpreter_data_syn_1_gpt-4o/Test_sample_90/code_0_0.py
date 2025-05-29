from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = -63*X**6 - 10*X**4 - 24*X**3

# Perform the integration
integral = integrate(function, X)

# Print the result
print(integral)