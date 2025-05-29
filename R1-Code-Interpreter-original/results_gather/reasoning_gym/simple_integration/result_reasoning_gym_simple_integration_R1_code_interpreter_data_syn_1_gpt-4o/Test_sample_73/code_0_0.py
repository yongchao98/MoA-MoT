from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -45*x**8 + 8*x**7 + 42*x**5 + 9*x**2/2

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)