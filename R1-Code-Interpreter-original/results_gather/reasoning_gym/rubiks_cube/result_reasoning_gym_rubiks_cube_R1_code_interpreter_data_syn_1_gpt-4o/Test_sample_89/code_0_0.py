# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'O', 'O', 'O', 'Y', 'Y', 'B', 'R', 'R', 'R',
    'B', 'R', 'W', 'G', 'G', 'G', 'Y', 'O', 'G',
    'W', 'W', 'W', 'Y', 'R', 'W', 'G', 'G', 'G',
    'Y', 'O', 'W', 'B', 'B', 'B', 'Y', 'Y', 'Y',
    'B', 'O', 'W', 'B', 'B', 'B', 'Y', 'R', 'G',
    'O', 'G', 'R', 'O', 'W', 'R', 'O', 'W', 'R'
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to Singmaster notation without combining steps
singmaster_solution = ' '.join(move.replace('2', ' ').replace('\'', ' \'') for move in solution)

# Print the solution
print(singmaster_solution)