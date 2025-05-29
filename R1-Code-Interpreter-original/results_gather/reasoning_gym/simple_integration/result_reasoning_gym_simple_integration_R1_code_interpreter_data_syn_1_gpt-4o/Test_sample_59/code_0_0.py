from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -56*x**7 + 20*x**3/9

# Perform the integration
integral = integrate(function, x)

# Print the result
print(integral)