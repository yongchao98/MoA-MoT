from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = (11*X**9)/7 + 30*X**4

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)