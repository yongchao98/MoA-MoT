# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the current state of the cube in a format that the solver can understand
# The state is represented as a string of 54 characters, each representing a color
# The order is: U (up), R (right), F (front), D (down), L (left), B (back)
# Colors: W (white), R (red), B (blue), O (orange), G (green), Y (yellow)

cube_state = (
    "YRWGYGYOW"
    "BWGOBRBYG"
    "RYORRGOGRB"
    "OWYBWBWOY"
)

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Convert the solution to the required format (no combined steps)
solution_steps = []
for move in solution:
    if len(move) == 2 and move[1] == '2':
        solution_steps.append(move[0])
        solution_steps.append(move[0])
    else:
        solution_steps.append(move)

# Print the solution
print(' '.join(solution_steps))