# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the current state of the cube in a string format
# The cube state is represented as a string of 54 characters, each representing a color
# The order is: UUUUUUUUURRRRRRRRRFFFFFFFFFDDDDDDDDDLLLLLLLLLBBBBBBBBB
# Using the given cube state:
cube_state = (
    "OOOYYYRRR"  # U face
    "YRWGGGYOW"  # R face
    "YRWGGGYOW"  # F face
    "GGGYOWBBY"  # D face
    "OWROWR"     # L face
    "BBBYRWR"    # B face
)

# Solve the cube using the library
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))