from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 9*x**8 - (32/5)*x**7 - 35*x**6

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)