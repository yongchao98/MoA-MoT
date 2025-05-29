from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = -20*X**(4/7) + 12*X

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)