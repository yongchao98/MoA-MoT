from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -6*x**5 - 24*x**3

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)