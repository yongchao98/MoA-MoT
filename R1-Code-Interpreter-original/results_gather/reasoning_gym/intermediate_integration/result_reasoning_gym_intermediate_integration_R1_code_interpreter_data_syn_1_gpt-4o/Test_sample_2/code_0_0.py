from sympy import symbols, integrate, exp

# Define the variable and the function
x = symbols('x')
function = -2*x*exp(x)

# Calculate the indefinite integral
integral_result = integrate(function, x)

# Print the result
print(integral_result)