from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = (9/2)*x**8 - (14/3)*x**5 - 20*x**4 + 10*x

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)