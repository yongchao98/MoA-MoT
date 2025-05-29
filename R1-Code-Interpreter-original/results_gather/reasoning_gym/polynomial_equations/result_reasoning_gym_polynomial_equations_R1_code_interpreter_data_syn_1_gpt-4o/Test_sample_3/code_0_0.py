import sympy as sp

# Define the variable
t = sp.symbols('t')

# Define the equation
equation = 7*t**2

# Solve the equation
solutions = sp.solve(equation, t)

# Format the solutions as decimal values rounded to 4 decimal places
formatted_solutions = [round(float(sol), 4) for sol in solutions]

# Print the formatted solutions
print(formatted_solutions)