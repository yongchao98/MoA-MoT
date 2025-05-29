from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'G', 'Y', 'Y', 'G', 'Y', 'Y', 'G', 'O', 'G',
    'Y', 'Y', 'Y', 'R', 'G', 'R', 'W', 'O', 'O', 'B', 'B', 'O',
    'R', 'R', 'B', 'W', 'G', 'R', 'W', 'O', 'O', 'B', 'B', 'Y',
    'R', 'R', 'B', 'W', 'G', 'R', 'W', 'G', 'G', 'O', 'O', 'Y',
    'O', 'R', 'B', 'B', 'W', 'W', 'B', 'W', 'W'
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to the required format (no combined steps)
formatted_solution = ' '.join(move.replace('2', ' ').replace('\'', ' \'') for move in solution)

print(formatted_solution)