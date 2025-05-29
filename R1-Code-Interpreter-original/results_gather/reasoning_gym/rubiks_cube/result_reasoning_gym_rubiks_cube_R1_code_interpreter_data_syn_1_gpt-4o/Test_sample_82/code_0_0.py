# Importing the necessary library
from rubik_solver import utils

# Define the initial state of the cube
# The state is represented as a string with each face's colors in order:
# U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The colors are represented by their initials: O, Y, B, R, G, W
cube_state = (
    "OOO"  # Up face
    "YYB"
    "BRB"
    "YRR"  # Right face
    "YGO"
    "WWG"
    "YRW"  # Front face
    "RGY"
    "OOG"
    "GBW"  # Down face
    "RGW"
    "RWG"
    "YRW"  # Left face
    "RGY"
    "OOG"
    "GBW"  # Back face
    "RGW"
    "RWG"
)

# Solve the cube using the Kociemba's algorithm
solution = utils.solve(cube_state)

# Print the solution in Singmaster notation
print(' '.join(solution))