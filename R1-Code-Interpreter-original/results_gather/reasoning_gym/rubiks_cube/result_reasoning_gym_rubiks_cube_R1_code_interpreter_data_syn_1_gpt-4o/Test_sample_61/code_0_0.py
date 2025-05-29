from rubik_solver import utils

# Define the cube state
cube_state = [
    'Y', 'W', 'R', 'O', 'Y', 'W', 'O', 'Y', 'Y',
    'R', 'Y', 'W', 'B', 'G', 'G', 'O', 'O', 'G',
    'Y', 'R', 'B', 'G', 'R', 'W', 'B', 'G', 'Y',
    'B', 'O', 'B', 'O', 'B', 'W', 'R', 'G', 'G',
    'Y', 'W', 'G', 'O', 'R', 'Y'
]

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Convert the solution to Singmaster notation without combining steps
singmaster_solution = ' '.join(move.replace('2', ' ').replace('\'', ' \'') for move in solution)

print(singmaster_solution)