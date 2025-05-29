from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = [
    'G', 'O', 'R', 'Y', 'Y', 'G', 'O', 'O', 'G',
    'Y', 'B', 'G', 'W', 'W', 'R', 'W', 'W', 'B', 'W', 'B', 'O',
    'Y', 'R', 'Y', 'R', 'G', 'B', 'W', 'O', 'O', 'G', 'B', 'O',
    'Y', 'R', 'Y', 'R', 'G', 'B', 'W', 'G', 'G', 'R', 'B', 'O',
    'B', 'R', 'O', 'W', 'W', 'Y', 'B', 'R', 'Y'
]

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))