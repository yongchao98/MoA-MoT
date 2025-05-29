from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = -54*X**(8/7) + 24*X**7 + 63/2*X**6 + 28*X**3

# Integrate the function
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)