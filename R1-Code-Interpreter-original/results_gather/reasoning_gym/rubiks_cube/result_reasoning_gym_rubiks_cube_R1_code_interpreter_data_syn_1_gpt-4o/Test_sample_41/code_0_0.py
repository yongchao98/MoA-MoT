# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'G', 'W', 'R', 'O', 'Y', 'R', 'O', 'G', 'R',
    'Y', 'B', 'B', 'Y', 'R', 'G', 'Y', 'Y', 'W',
    'G', 'B', 'O', 'Y', 'R', 'R', 'W', 'G', 'G',
    'O', 'O', 'Y', 'B', 'B', 'G', 'B', 'R', 'R',
    'W', 'G', 'G', 'O', 'O', 'B', 'W', 'O', 'R',
    'B', 'W', 'W', 'B', 'W', 'W', 'Y', 'Y', 'O'
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to a string with spaces between moves
solution_str = ' '.join(solution)

# Print the solution
print(solution_str)