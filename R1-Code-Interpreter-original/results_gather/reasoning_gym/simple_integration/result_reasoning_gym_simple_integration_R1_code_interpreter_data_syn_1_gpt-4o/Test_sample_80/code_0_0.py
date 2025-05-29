from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
function = 35*X**9 - 35*X**6/4 + 67*X**5/12 - 30*X**2

# Integrate the function
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)