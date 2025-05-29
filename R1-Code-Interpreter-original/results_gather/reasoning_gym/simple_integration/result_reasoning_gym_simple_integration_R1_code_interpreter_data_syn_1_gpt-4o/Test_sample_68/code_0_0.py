from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = -21*X**(5/2) - (2*X**3)/3 + 18*X

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)