from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 24*x**8 - 30*x**5 - 6

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)