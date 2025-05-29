# Importing a Rubik's cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = [
    'R', 'Y', 'Y', 'G', 'Y', 'O', 'Y', 'Y', 'G',
    'Y', 'Y', 'O', 'B', 'B', 'R', 'W', 'B', 'G', 'O', 'O', 'G',
    'G', 'R', 'R', 'W', 'G', 'B', 'R', 'O', 'O', 'W', 'B', 'W',
    'B', 'R', 'R', 'B', 'W', 'R', 'W', 'O', 'O', 'W', 'Y', 'O',
    'Y', 'B', 'B', 'G', 'W', 'G', 'W', 'R', 'G'
]

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to the required format
formatted_solution = ' '.join(solution)

# Print the solution
print(formatted_solution)