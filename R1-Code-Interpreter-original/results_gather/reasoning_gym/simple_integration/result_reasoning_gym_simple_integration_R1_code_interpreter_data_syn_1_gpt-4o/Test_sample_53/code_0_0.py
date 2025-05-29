from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function
f = 63*x**8 - (32/5)*x**7 + 1/4

# Calculate the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)