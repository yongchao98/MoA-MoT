from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = -8*X**7 + 146*X**6/9

# Integrate the function
integral_result = integrate(function, X)

# Print the result
print(integral_result)