from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function to integrate
function = 48*X**5 - 10/9

# Compute the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)