from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = -5*X**9 + 63*X**8/4 - 12*X**5 - 15*X**2

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)