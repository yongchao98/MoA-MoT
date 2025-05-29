from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
f = 3 * (6*x - 9)**3

# Calculate the indefinite integral
indefinite_integral = integrate(f, x)

# Print the result
print(indefinite_integral)