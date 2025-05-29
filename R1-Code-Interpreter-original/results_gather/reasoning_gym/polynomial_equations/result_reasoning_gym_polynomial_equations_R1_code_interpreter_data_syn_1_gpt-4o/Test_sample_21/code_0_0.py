from sympy import symbols, Eq, solve

# Define the variable
v = symbols('v')

# Define the equation
equation = Eq(-85 * v, 0)

# Solve the equation
solution = solve(equation, v)

# Print the solution
print(solution)