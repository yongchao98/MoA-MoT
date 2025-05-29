from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'O', 'R', 'W', 'B', 'Y', 'B', 'B', 'R', 'B',
    'G', 'W', 'W', 'R', 'G', 'R', 'Y', 'Y', 'B',
    'O', 'B', 'Y', 'R', 'R', 'B', 'O', 'G', 'Y',
    'G', 'O', 'O', 'W', 'B', 'W', 'W', 'W', 'W',
    'R', 'G', 'B', 'O', 'Y', 'Y', 'G', 'O', 'O',
    'G', 'O', 'Y', 'G', 'W', 'R', 'G', 'Y', 'R'
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Print the solution
print(' '.join(solution))