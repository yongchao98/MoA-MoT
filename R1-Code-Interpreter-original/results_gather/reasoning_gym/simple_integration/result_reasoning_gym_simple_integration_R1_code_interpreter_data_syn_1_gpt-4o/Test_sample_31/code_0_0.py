from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = 55*X**9 + 18*X**5 - 8*X**3 - 8*X/3

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)