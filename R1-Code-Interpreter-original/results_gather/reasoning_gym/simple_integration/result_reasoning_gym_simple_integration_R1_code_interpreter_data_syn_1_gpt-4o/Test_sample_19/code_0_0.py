from sympy import symbols, integrate

# Define the variable
x = symbols('x')

# Define the expressions
expr1 = -216 * x**(8/5)
expr2 = 35 * x**(3/4)

# Integrate the expressions
integral1 = integrate(expr1, x)
integral2 = integrate(expr2, x)

# Print the results
print(integral1)
print(integral2)