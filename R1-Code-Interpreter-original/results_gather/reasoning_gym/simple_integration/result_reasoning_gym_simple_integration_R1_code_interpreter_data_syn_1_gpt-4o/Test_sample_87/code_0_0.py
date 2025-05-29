from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the function
function = (14*x**6)/3 + (15*x**4)/2 + 24*x**3

# Calculate the antiderivative
antiderivative = integrate(function, x)

# Print the result
print(antiderivative)