from rubik_solver import utils

# Define the current state of the cube using the color notation
# The cube is represented as a string of 54 characters, each representing a color
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# Colors: W (White), Y (Yellow), R (Red), O (Orange), G (Green), B (Blue)

cube_state = (
    "GWW"
    "OYO"
    "GOO"
    "YBYRGGW"
    "WGRRO"
    "YR"
    "YRGBWO"
    "RBB"
    "O"
    "BYRBGO"
    "YGBR"
    "YO"
    "YWB"
    "GWR"
    "WBW"
)

# Solve the cube using the library
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))