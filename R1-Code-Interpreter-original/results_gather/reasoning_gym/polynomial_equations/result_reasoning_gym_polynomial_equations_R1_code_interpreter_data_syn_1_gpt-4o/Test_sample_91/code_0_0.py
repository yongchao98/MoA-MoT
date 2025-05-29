# Solving the equation -35*b**2 - 5*b = 0
import sympy as sp

# Define the variable
b = sp.symbols('b')

# Define the equation
equation = -35*b**2 - 5*b

# Solve the equation
solutions = sp.solve(equation, b)

# Convert solutions to decimal and round to 4 decimal places
decimal_solutions = [round(float(sol), 4) for sol in solutions]

# Print the solutions
print(decimal_solutions)