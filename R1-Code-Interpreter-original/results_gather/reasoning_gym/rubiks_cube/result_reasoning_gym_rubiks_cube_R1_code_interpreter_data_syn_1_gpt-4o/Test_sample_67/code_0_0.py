# Importing a library to solve the Rubik's Cube
from rubik_solver import utils

# Define the scrambled state of the cube
# The state is represented as a string of 54 characters, each representing a color
# The order is: UUUUUUUUURRRRRRRRRFFFFFFFFFDDDDDDDDDLLLLLLLLLBBBBBBBBB
# Using the given scrambled state:
scrambled_state = (
    "BOYYYGGRB"  # U face
    "YGOWGORYOR"  # R face
    "YRYBGOYOBR"  # F face
    "GGBRWWGWOY"  # D face
    "WOOOWWBBBR"  # L face
    "GRBRWWBBBR"  # B face
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))