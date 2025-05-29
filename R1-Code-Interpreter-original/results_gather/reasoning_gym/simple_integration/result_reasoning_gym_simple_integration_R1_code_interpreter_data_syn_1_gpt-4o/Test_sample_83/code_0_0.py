from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 80*x**9 + 8*x**7 - 18*x - 3

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)