from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -9*x**8/2 - 21*x**2 + 4*x + 1

# Compute the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)