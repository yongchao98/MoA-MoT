# Solving the equation -164*g = 0
from sympy import symbols, Eq, solve

# Define the variable
g = symbols('g')

# Define the equation
equation = Eq(-164*g, 0)

# Solve the equation
solution = solve(equation, g)

# Print the solution
print(solution)