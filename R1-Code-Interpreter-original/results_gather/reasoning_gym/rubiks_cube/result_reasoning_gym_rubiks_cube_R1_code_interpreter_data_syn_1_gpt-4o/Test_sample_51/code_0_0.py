from rubik_solver import utils

# Define the cube configuration in the required format
# The cube is represented as a string of 54 characters, each representing a color
# The order is: UUUUUUUUURRRRRRRRRFFFFFFFFFDDDDDDDDDLLLLLLLLLBBBBBBBBB
# Using the given configuration:
# U: R G Y G Y Y G B O
# R: G Y Y R O G W R G
# F: O W W R R B W G O
# D: W R G O W W O G B
# L: R R B W G O Y O Y
# B: B B B R O Y R O Y

cube_configuration = (
    "RGYGYYGBO"  # U
    "GYYROGWRG"  # R
    "OWWRRBWGO"  # F
    "WRGOWWOGB"  # D
    "RRBWGOYOY"  # L
    "BBBROYROY"  # B
)

# Solve the cube
solution = utils.solve(cube_configuration, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))