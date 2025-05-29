from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -x**7 - 32*x**3 + 5

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)