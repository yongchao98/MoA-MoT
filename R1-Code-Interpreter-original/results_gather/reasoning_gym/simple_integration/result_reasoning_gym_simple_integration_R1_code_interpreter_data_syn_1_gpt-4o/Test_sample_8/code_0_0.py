from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function
function = -70*x**9 + 30*x**5 + 20*x**3 + 4/3

# Calculate the indefinite integral
indefinite_integral = integrate(function, x)

# Print the result
print(indefinite_integral)