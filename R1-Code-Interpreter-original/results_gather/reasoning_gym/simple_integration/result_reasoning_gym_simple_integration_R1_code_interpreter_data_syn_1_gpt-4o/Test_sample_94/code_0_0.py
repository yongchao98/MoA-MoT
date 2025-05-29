from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = -24*X**7 - 9*X**2

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)