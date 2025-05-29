from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = 30*X**4 - 7*X**3 + 46*X**2/3

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)