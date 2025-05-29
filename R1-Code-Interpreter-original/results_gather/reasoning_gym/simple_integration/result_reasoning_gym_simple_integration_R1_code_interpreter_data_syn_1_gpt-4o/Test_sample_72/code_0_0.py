from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 42*x**5 + 3*x**2

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)