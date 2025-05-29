# Importing the necessary library
from rubik_solver import utils

# Define the current state of the cube
# The state is represented as a string with each face's colors in order:
# U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The colors are represented by their initials: W (White), Y (Yellow), R (Red), O (Orange), G (Green), B (Blue)

cube_state = (
    "OYO"  # Up face
    "WYB"
    "OOO"
    "BRB"  # Right face
    "YGY"
    "BRY"
    "YGY"  # Front face
    "GOR"
    "YOR"
    "RGW"  # Down face
    "RWW"
    "RYR"
    "BRB"  # Left face
    "YGY"
    "BRY"
    "GOG"  # Back face
    "WBW"
    "WBW"
)

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))