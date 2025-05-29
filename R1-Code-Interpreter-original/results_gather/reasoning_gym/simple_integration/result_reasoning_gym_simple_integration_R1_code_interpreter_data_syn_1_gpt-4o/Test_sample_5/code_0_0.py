from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 7*x**7 - 2*x

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)