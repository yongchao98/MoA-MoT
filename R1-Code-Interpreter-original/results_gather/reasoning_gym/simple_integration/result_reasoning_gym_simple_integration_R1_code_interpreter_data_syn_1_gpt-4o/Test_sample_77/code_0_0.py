from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -45/4*x**8 - 104*x**7 + 1

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)