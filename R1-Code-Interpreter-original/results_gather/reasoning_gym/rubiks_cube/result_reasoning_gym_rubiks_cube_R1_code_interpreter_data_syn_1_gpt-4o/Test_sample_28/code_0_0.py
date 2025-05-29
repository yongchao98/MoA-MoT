# Importing the necessary library
from rubik_solver import utils

# Define the current state of the cube using the given configuration
# The cube is represented in a string format where each face is represented by a single character:
# U (Up), D (Down), F (Front), B (Back), L (Left), R (Right)
# The colors are represented as: W (White), Y (Yellow), R (Red), G (Green), B (Blue), O (Orange)

# The cube state is given in the following order: U, R, F, D, L, B
# We need to convert the given configuration into this format

cube_state = (
    "WWWWWWWWW"  # U face (White)
    "RRRRRRRRR"  # R face (Red)
    "YYYYYYYYY"  # F face (Yellow)
    "YYYYYYYYY"  # D face (Yellow)
    "GGGGGGGGG"  # L face (Green)
    "BBBBBBBBB"  # B face (Blue)
)

# Solve the cube using the library
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))