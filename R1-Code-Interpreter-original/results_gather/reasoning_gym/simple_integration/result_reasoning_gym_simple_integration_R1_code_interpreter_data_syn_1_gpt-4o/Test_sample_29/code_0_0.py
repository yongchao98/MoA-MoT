from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -56*x**7 - 48*x**5 - 30*x**4

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)