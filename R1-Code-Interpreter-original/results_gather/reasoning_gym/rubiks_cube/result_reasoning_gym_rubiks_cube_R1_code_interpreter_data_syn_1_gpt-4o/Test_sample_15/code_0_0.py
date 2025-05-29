from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = [
    'R', 'Y', 'G', 'B', 'Y', 'G', 'G', 'W', 'W',
    'W', 'O', 'O', 'Y', 'B', 'B', 'O', 'R', 'R', 'W', 'G', 'B',
    'G', 'R', 'O', 'Y', 'G', 'G', 'O', 'O', 'R', 'Y', 'B', 'W',
    'B', 'O', 'R', 'Y', 'B', 'B', 'O', 'R', 'W', 'G', 'R', 'R',
    'G', 'Y', 'Y', 'W', 'W', 'B', 'Y', 'W', 'O'
]

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to the required format
formatted_solution = ' '.join(solution).replace('2', ' 2').replace("'", " '")

print(formatted_solution)