from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = 30*X**5 - 8*X**3 - 10*X

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)