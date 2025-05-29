from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 70*x**6 + 12*x**2/5

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)