from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 90*x**8 - 7*x**6 + 30*x**5 - 35*x**4

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)