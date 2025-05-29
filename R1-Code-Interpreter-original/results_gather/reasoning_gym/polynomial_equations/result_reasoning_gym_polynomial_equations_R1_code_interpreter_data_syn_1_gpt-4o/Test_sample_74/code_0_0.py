from sympy import symbols, Eq, solve

# Define the variable
y = symbols('y')

# Define the equation
equation = Eq(76 * y**2, 0)

# Solve the equation
solutions = solve(equation, y)

# Print the solutions
print(solutions)