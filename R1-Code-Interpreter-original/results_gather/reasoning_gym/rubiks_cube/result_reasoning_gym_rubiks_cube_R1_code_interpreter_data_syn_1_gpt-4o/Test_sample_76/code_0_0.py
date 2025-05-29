# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
# The cube is represented in a string format where each face is listed in order:
# U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# Each face is represented by 9 characters (3x3), and the colors are represented by their initials:
# W (White), Y (Yellow), G (Green), B (Blue), R (Red), O (Orange)

scrambled_state = (
    "WYY"  # Up face
    "GYY"
    "GYY"
    "GYR"  # Right face
    "ROO"
    "BBB"
    "RRO"  # Front face
    "WGG"
    "WGG"
    "WOO"  # Down face
    "BRR"
    "BWW"
    "YWW"  # Left face
    "GRR"
    "GRR"
    "WOO"  # Back face
    "BOO"
    "BBO"
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))