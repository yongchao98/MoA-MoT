from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -56*x**6 + 15*x**5 + 15*x**4/2

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)