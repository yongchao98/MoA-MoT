from rubik_solver import utils

# Define the scrambled cube state
scrambled_state = (
    'G', 'O', 'Y', 'R', 'Y', 'B', 'R', 'Y', 'O',
    'O', 'B', 'B', 'W', 'R', 'Y', 'G', 'Y', 'B', 'R', 'B', 'W',
    'O', 'R', 'R', 'G', 'G', 'Y', 'G', 'O', 'W', 'R', 'B', 'W',
    'Y', 'G', 'Y', 'G', 'O', 'O', 'B', 'B', 'W', 'R', 'O', 'O',
    'R', 'G', 'W', 'W', 'W', 'W', 'B', 'Y', 'G'
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to a string with each move separated by a space
solution_str = ' '.join(solution)

# Print the solution
print(solution_str)