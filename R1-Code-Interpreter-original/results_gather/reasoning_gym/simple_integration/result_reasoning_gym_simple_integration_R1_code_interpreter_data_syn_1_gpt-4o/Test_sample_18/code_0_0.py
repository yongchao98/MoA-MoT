from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = -20*X**9 + 48*X**5/5 + 18*X/5

# Perform the integration
integral = integrate(function, X)

# Print the result
print(integral)