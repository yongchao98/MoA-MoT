# Importing a library to solve the Rubik's Cube
from rubik_solver import utils

# Define the cube state in a format that the solver can understand
# The cube state is given in a specific order: U, R, F, D, L, B
# Each face is represented by a string of 9 characters, where each character is a color
# The colors are represented by their initials: W (white), R (red), B (blue), O (orange), G (green), Y (yellow)

# Given cube state:
# U: G G G
#    O Y Y
#    O Y Y
# L: Y B B
#    Y R R
#    W W W
# F: O O O
#    Y B B
#    B R R
# R: G G W
#    O O O
#    Y R R
# B: G G R
#    B R R
#    B W W
# D: B W W
#    B W W
#    Y R R

# Convert the given cube state into a string format
cube_state = (
    "GGGOYYOYY"  # U
    "YBBYRRWWW"  # L
    "OOOYBBBRR"  # F
    "GGWOOOYRR"  # R
    "GGRBRRWWB"  # B
    "BWWBWWYRR"  # D
)

# Solve the cube using the library
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))