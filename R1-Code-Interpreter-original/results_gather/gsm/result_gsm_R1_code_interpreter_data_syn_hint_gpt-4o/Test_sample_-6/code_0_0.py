from sympy import symbols, Eq, solve

# Define the variable
r = symbols('r')

# Define the equation based on the problem statement
equation = Eq(3*r + 1801110*(r + 2), 42)

# Solve the equation
red_chickens = solve(equation, r)

# Print the result
print(red_chickens)