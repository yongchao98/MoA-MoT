from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = [
    'R', 'R', 'W',
    'Y', 'Y', 'Y',
    'W', 'W', 'W',
    'Y', 'B', 'G', 'R', 'R', 'R', 'B', 'G', 'B', 'O', 'G', 'G',
    'W', 'R', 'R', 'B', 'G', 'Y', 'O', 'O', 'B', 'O', 'B', 'B',
    'G', 'R', 'R', 'B', 'G', 'Y', 'O', 'O', 'B', 'O', 'W', 'W',
    'Y', 'W', 'G',
    'Y', 'W', 'G',
    'O', 'O', 'Y'
]

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to a string with each move separated by a space
solution_str = ' '.join(solution)

# Print the solution
print(solution_str)