import sympy as sp

# Define the variable
k = sp.symbols('k')

# Define the equation
equation = -k**2 + 148*k

# Solve the equation
solutions = sp.solve(equation, k)

# Format the solutions to four decimal places
formatted_solutions = [f"{sol.evalf():.4f}" for sol in solutions]

# Print the formatted solutions
print(", ".join(formatted_solutions))