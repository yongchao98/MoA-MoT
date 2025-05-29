from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 50*x**4 - 2*x**3/5 + 6*x**2 - 32*x/5

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)