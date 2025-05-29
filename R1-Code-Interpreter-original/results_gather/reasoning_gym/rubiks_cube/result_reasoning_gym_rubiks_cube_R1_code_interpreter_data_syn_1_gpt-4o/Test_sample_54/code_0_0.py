from rubik_solver import utils

# Define the current state of the cube using the color notation
# The cube state is represented as a string of 54 characters, each representing a color
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The given cube state is:
# U: G R R G Y Y G R R
# R: W W W R R R W W W
# F: O G G O G G O G G
# D: Y O Y Y O Y Y O Y
# L: B B R B B Y B B R
# B: B O O B W W B O O

cube_state = (
    "GRRGYYGRR"  # U
    "WWW"        # R
    "RRR"        # R
    "WWW"        # R
    "OGGOGGOGG"  # F
    "YOYYOYYOY"  # D
    "BBRBBYBBR"  # L
    "BOOBWWBOO"  # B
)

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))