import sympy as sp

# Define the variable and equation
k = sp.symbols('k')
equation = -8*k**3 - 47*k

# Solve the equation
solutions = sp.solve(equation, k)

# Convert solutions to decimal and round to 4 decimal places
decimal_solutions = [float(sol.evalf()) for sol in solutions]
rounded_solutions = [round(sol, 4) for sol in decimal_solutions]

# Print the solutions
print(rounded_solutions)