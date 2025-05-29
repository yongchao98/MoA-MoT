from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = 100*X**9 + 90*X**8

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)