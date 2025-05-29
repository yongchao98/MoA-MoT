from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function to integrate
function = -8*X**(5/3) + 38*X**3

# Perform the integration
integral = integrate(function, X)

# Print the result
print(integral)