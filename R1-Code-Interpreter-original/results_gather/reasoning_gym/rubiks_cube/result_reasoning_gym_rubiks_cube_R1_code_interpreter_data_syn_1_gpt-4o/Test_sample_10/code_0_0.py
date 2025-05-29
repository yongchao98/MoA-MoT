from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'G', 'Y', 'W', 'G', 'Y', 'W', 'R', 'R', 'R',
    'R', 'R', 'B', 'W', 'W', 'W', 'G', 'O', 'O', 'G', 'B', 'Y',
    'R', 'R', 'W', 'G', 'G', 'G', 'Y', 'O', 'O', 'G', 'B', 'Y',
    'R', 'R', 'Y', 'B', 'B', 'B', 'W', 'O', 'O', 'G', 'B', 'Y',
    'O', 'O', 'O', 'B', 'W', 'Y', 'B', 'W', 'Y'
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to the required format
formatted_solution = ' '.join(solution).replace('2', ' ').replace('\'', ' \'').replace('  ', ' ')

print(formatted_solution)