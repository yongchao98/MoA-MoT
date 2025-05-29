import sympy as sp

# Define the variable
v = sp.symbols('v')

# Define the equation
equation = 16*v**3 + 24*v**2

# Solve the equation
solutions = sp.solve(equation, v)

# Convert solutions to decimal and round to 4 decimal places
decimal_solutions = [float(sol.evalf()) for sol in solutions]
rounded_solutions = [round(sol, 4) for sol in decimal_solutions]

# Print the solutions
print(rounded_solutions)