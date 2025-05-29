from rubik_solver import utils

# Define the cube configuration
# The cube is represented in a string format where each face is represented by a single letter
# The order of faces is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The given configuration is translated into this format
cube_configuration = (
    "RYG"  # U face
    "BYR"
    "WYR"
    "GRG"  # R face
    "GGO"
    "YBW"
    "GRO"  # F face
    "RGO"
    "RBB"
    "WWO"  # D face
    "WWO"
    "BBY"
    "GYO"  # L face
    "RGY"
    "YBB"
    "YBO"  # B face
    "YOW"
    "WGO"
)

# Solve the cube
solution = utils.solve(cube_configuration, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))