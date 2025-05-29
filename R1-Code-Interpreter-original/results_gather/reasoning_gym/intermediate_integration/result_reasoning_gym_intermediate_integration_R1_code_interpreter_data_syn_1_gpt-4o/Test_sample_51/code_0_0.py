from sympy import symbols, sin, cos, integrate

# Define the variable and the function
x = symbols('x')
u = 7*x + 1
f = 7*sin(u)**2*cos(u)

# Perform the integration
antiderivative = integrate(f, x)

# Output the result
print(antiderivative)