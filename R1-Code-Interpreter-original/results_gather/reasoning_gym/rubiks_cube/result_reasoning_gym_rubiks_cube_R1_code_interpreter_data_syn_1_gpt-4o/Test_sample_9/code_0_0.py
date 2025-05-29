from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'B', 'Y', 'B', 'B', 'Y', 'B', 'Y', 'Y', 'G',
    'O', 'O', 'G', 'O', 'O', 'R', 'W', 'R', 'R', 'Y', 'G', 'Y',
    'G', 'R', 'R', 'G', 'G', 'W', 'O', 'O', 'O', 'G', 'B', 'W',
    'G', 'B', 'B', 'R', 'R', 'G', 'R', 'R', 'W', 'O', 'B', 'W',
    'W', 'W', 'Y', 'W', 'W', 'Y', 'O', 'Y', 'B'
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to the required format
formatted_solution = ' '.join(solution).replace('2', ' 2').replace('\'', ' \'')
print(formatted_solution)