from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function to integrate
function = 10*X**4 + 14*X**3/5

# Perform the integration
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)