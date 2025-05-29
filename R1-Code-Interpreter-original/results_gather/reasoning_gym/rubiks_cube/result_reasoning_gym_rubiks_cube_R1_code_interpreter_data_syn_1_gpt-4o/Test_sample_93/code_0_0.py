# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Defining the cube configuration
# The configuration is given in the order: U, R, F, D, L, B
# Each face is represented by a string of 9 characters, where:
# W = White, Y = Yellow, R = Red, O = Orange, G = Green, B = Blue

cube_configuration = (
    "WOO"  # U face
    "BYB"
    "WYB"
    "BRR"  # R face
    "YGY"
    "YRR"
    "BYY"  # F face
    "OGY"
    "OOG"
    "GWW"  # D face
    "GWW"
    "YGG"
    "BYB"  # L face
    "BRR"
    "BRR"
    "OOG"  # B face
    "WWR"
    "BBR"
)

# Solving the cube
solution = utils.solve(cube_configuration, 'Kociemba')

# Printing the solution in Singmaster notation
print(' '.join(solution))