# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled cube state
scrambled_state = (
    'W', 'G', 'W', 'O', 'Y', 'Y', 'O', 'B', 'B',
    'O', 'Y', 'Y', 'G', 'Y', 'Y', 'O', 'G', 'R',
    'B', 'O', 'B', 'O', 'R', 'R', 'Y', 'G', 'R',
    'B', 'O', 'W', 'R', 'B', 'B', 'Y', 'O', 'R',
    'B', 'W', 'O', 'W', 'W', 'W', 'R', 'R', 'R',
    'Y', 'B', 'G', 'W', 'W', 'G', 'G', 'G', 'G'
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to Singmaster notation without combining steps
singmaster_solution = ' '.join(move.replace('2', ' ').replace('\'', ' \'') for move in solution).replace('  ', ' ')

# Print the solution
print(singmaster_solution)