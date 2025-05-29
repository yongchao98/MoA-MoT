from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = 24*x**7 - 28*x**6/5 + 33*x**5

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)