# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the current state of the Rubik's Cube
# The state is represented as a string of 54 characters, each representing a color
# The order is: UUUUUUUUURRRRRRRRRFFFFFFFFFDDDDDDDDDLLLLLLLLLBBBBBBBBB
# Using the provided cube state:
cube_state = (
    "ROYGBYWRO"
    "YOBYGRWRY"
    "RGBOWBYBW"
    "YORWWWORG"
)

# Solve the cube using the solver
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))