# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the current state of the cube using the color notation
# The cube is represented as a string of 54 characters, each representing a sticker
# The order is: UUUUUUUUURRRRRRRRRFFFFFFFFFDDDDDDDDDLLLLLLLLLBBBBBBBBB
# Using the given configuration:
cube_state = (
    "GGBYYOYRO"
    "YRRGWGRYO"
    "YRRGGBBOW"
    "WBYYGWBWR"
    "OBBRRWGWW"
    "OOGWBBRRB"
)

# Solve the cube using the library
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))