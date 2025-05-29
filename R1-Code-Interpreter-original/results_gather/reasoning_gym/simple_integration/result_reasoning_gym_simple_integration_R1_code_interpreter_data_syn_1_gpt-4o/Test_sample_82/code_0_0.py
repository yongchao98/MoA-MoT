from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -16*x**7 - 18*x**5 + 2

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)