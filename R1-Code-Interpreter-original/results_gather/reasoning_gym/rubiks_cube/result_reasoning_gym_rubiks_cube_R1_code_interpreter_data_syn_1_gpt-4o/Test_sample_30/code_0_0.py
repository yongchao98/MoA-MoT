from rubik_solver import utils

# Define the cube state in a format that the solver can understand
# The cube state is given in the order: U, R, F, D, L, B
# Each face is represented as a string of 9 characters, where each character is the color of a sticker
# The colors are represented by their initials: W (white), R (red), B (blue), O (orange), G (green), Y (yellow)

cube_state = (
    "YYRYOYOWRRYO"  # U face
    "RGBRGBYOWOWG"  # R face
    "BBBRGBOWGBWG"  # F face
    "WRYGYWGWOW"    # D face
    "RRBRYOGGWGB"   # L face
    "BOOYYOOGBBW"   # B face
)

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))