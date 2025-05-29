from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = -21*X**6 - 138*X**5/5 - 9/4

# Integrate the function
integral = integrate(function, X)

# Print the result
print(integral)