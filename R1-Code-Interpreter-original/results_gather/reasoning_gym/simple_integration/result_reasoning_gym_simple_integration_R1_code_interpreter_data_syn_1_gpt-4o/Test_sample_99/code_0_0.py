from sympy import symbols, integrate

# Define the variable and the function
X = symbols('X')
function = 20*X**9 - 144*X**8/5 + 35*X**4

# Calculate the indefinite integral
indefinite_integral = integrate(function, X)

# Print the result
print(indefinite_integral)