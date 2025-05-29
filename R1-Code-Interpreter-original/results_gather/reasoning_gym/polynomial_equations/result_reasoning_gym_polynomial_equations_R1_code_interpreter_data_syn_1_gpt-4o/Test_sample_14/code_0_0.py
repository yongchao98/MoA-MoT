import sympy as sp

# Define the variable
c = sp.symbols('c')

# Define the equation
equation = 106 * c**2

# Solve the equation
solutions = sp.solve(equation, c)

# Print the solutions rounded to 4 decimal places
print([round(float(sol), 4) for sol in solutions])