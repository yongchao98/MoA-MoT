from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -8*x**9 + 49*x**6/10 - 3*x**5/5 - 1

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)