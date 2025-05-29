from sympy import symbols, integrate

# Define the variable and the function
x = symbols('x')
function = -9*x**8 - 56*x**7 + 6*x**5 + 16*x**3 + 6

# Perform the integration
integral = integrate(function, x)

# Print the result
print(integral)