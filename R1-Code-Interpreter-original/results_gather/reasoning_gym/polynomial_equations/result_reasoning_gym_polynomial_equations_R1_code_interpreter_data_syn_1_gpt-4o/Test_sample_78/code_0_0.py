from sympy import symbols, solve

# Define the variable
o = symbols('o')

# Define the equation
equation = 144 * o**3

# Solve the equation
solutions = solve(equation, o)

# Print the solutions
print(solutions)