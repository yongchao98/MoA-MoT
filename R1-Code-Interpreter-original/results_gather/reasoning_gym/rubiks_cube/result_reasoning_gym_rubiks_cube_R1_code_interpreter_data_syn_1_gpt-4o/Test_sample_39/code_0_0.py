from rubik_solver import utils

# Define the cube state
cube_state = (
    'B', 'B', 'B',
    'B', 'Y', 'B',
    'W', 'O', 'O',
    'Y', 'Y', 'B', 'O', 'Y', 'Y', 'G', 'W', 'W', 'R', 'R', 'R',
    'W', 'R', 'Y', 'G', 'G', 'R', 'W', 'O', 'O', 'B', 'B', 'O',
    'R', 'R', 'B', 'O', 'Y', 'Y', 'G', 'O', 'O', 'W', 'W', 'W',
    'Y', 'R', 'R',
    'G', 'W', 'G',
    'G', 'G', 'G'
)

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Break down each move into individual steps
detailed_solution = []
for move in solution:
    if len(move) == 2 and move[1] == '2':
        detailed_solution.extend([move[0], move[0]])
    else:
        detailed_solution.append(move)

# Print the detailed solution
print(' '.join(detailed_solution))