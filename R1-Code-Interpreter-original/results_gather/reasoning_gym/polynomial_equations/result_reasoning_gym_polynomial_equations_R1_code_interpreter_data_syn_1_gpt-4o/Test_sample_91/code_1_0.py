# Solving the equation -35*b**2 - 5*b = 0
import sympy as sp

# Define the variable
b = sp.symbols('b')

# Define the equation
equation = -35*b**2 - 5*b

# Solve the equation
solutions = sp.solve(equation, b)

# Custom rounding function to round down when the 5th decimal is 5 or greater
def custom_round(value, decimals):
    factor = 10 ** decimals
    return (int(value * factor) / factor)

# Convert solutions to decimal and apply custom rounding to 4 decimal places
decimal_solutions = [custom_round(float(sol), 4) for sol in solutions]

# Print the solutions
print(decimal_solutions)