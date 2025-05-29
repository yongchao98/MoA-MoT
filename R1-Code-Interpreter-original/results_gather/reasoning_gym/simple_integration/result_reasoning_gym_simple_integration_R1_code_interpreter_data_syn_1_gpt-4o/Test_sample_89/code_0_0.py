from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = (45/2)*x**8 + (7/2)*x**6

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)