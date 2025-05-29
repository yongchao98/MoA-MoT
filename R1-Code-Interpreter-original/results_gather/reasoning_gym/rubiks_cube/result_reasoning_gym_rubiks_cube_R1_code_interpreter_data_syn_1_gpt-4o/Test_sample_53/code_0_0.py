# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = [
    'G', 'G', 'R',  # Top face
    'Y', 'Y', 'B',
    'O', 'O', 'R',
    'Y', 'R', 'B', 'Y', 'Y', 'B', 'W', 'O', 'G', 'Y', 'O', 'O',  # Middle layer
    'Y', 'R', 'R', 'G', 'G', 'B', 'W', 'O', 'G', 'Y', 'B', 'B',
    'Y', 'O', 'R', 'G', 'G', 'B', 'W', 'R', 'W', 'G', 'B', 'B',
    'W', 'W', 'O',  # Bottom face
    'W', 'W', 'W',
    'R', 'R', 'O'
]

# Convert the scrambled state to a string format required by the solver
scrambled_string = ''.join(scrambled_state)

# Solve the cube
solution = utils.solve(scrambled_string, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))