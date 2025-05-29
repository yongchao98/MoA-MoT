from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = -16*X**7 + 12*X**5 - 4*X**3

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)