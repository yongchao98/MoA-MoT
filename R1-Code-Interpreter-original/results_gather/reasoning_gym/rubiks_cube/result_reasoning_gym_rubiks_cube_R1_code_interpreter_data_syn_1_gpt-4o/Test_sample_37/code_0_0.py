# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
# The state is represented in a string format where each face is represented by a single letter
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
scrambled_state = 'YYYYYYYYYRRRGGGOOOBBBRRRGGGOOOBBBGGGOOOBBBRRRWWWWWWWWW'

# Solve the cube using the solver
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to a string with each move separated by a space
solution_str = ' '.join(solution)

# Print the solution
print(solution_str)