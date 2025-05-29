# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = [
    'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'B', 'R', 'W',
    'G', 'G', 'W', 'R', 'G', 'O', 'G', 'B', 'B', 'R', 'R', 'R',
    'R', 'R', 'W', 'R', 'G', 'O', 'G', 'O', 'O', 'Y', 'B', 'B',
    'Y', 'B', 'W', 'R', 'W', 'O', 'G', 'G', 'W', 'O', 'O', 'O',
    'G', 'O', 'Y', 'W', 'W', 'W', 'B', 'B', 'B'
]

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to Singmaster notation without combining steps
singmaster_solution = ' '.join(move.replace('2', ' ').replace('\'', ' \'') for move in solution)

# Print the solution
print(singmaster_solution)