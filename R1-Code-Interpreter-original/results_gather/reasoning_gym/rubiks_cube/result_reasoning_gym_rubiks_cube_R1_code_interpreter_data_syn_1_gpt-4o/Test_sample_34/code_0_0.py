# Importing the necessary library
from rubik_solver import utils

# Defining the cube state in a format that the library can understand
# The cube state is given in a specific order: U, R, F, D, L, B
# We need to convert the given state into this format
cube_state = (
    "BBR"  # U
    "YYB"
    "YYB"
    "RRRGGYOO"  # R
    "YRRGGYOO"
    "YOYOOBWR"
    "GGOWWW"  # D
    "RRO"
    "GGY"  # L
    "WWB"
    "WWB"
    "OOW"  # B
    "GWB"
    "GWB"
)

# Solving the cube
solution = utils.solve(cube_state, 'Kociemba')

# Printing the solution in the required format
print(' '.join(solution))