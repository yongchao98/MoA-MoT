from sympy import symbols, exp, integrate

# Define the variable and function
x = symbols('x')
f = 3*x**3*exp(x)

# Compute the antiderivative
antiderivative = integrate(f, x)

# Print the result
print(antiderivative)