from sympy import symbols, integrate, log

# Define the variable
x = symbols('x')

# Define the function to integrate
f = -log(x)

# Perform the integration
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)