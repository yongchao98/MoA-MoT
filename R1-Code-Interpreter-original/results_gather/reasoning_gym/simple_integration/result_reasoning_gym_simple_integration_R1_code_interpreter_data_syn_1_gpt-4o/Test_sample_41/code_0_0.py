from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function
f = (81/8)*x**8 - 42*x**5 + (25/3)*x**4

# Calculate the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)