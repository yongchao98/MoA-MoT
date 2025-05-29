from sympy import symbols, integrate

# Define the variable
X = symbols('X')

# Define the function
f = 54*X**8 - (12*X**5)/5 - 25*X**4 - 1

# Calculate the indefinite integral
indefinite_integral = integrate(f, X)

# Print the result
print(indefinite_integral)