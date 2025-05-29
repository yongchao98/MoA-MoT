from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
f = 70*X**9 + 57*X**2/2 - 39*X/5

# Compute the indefinite integral
indefinite_integral = integrate(f, X)

# Print the result
print(indefinite_integral)