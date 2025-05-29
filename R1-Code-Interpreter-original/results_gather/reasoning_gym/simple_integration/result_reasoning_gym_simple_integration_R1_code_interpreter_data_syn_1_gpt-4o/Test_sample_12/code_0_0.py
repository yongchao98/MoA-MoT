from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function
function = -27*x**8 + 9/2*x**5 - 28*x**3 + 5*x**2 + 8*x

# Integrate the function
integral = integrate(function, x)

# Print the result
print(integral)