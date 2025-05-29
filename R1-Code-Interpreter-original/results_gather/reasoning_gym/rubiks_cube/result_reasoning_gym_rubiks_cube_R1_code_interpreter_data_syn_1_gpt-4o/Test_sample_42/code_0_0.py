from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = [
    'G', 'Y', 'G', 'G', 'Y', 'B', 'B', 'Y', 'B',  # Top face
    'R', 'R', 'O', 'W', 'G', 'Y', 'O', 'O', 'R', 'W', 'B', 'Y',  # Middle layer
    'R', 'R', 'O', 'W', 'G', 'Y', 'O', 'O', 'R', 'W', 'B', 'Y',  # Middle layer
    'R', 'R', 'O', 'W', 'B', 'Y', 'O', 'O', 'R', 'W', 'G', 'Y',  # Middle layer
    'G', 'W', 'G', 'B', 'W', 'G', 'B', 'W', 'B'   # Bottom face
]

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))