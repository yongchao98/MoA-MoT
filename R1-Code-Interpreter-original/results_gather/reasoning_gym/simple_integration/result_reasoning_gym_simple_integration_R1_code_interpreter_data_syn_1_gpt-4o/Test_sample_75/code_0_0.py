from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = 6*X**9 + 7*X**6/2 - 5/3

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)