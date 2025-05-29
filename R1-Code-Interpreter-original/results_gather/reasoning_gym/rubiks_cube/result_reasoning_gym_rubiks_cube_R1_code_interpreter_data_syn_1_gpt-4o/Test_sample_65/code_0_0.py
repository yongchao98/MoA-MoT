# Importing a library to solve the Rubik's Cube
from rubik_solver import utils

# Define the scrambled state of the cube
# The cube is represented as a string with each face's colors in order:
# U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
scrambled_state = (
    "ROOGYYYYY"  # U face
    "BRGORRGGW"  # R face
    "BRWRGGOWW"  # F face
    "GOWBWWWOO"  # D face
    "BRWRYBRGG"  # L face
    "WBBYYBBYY"  # B face
)

# Solve the cube using the Kociemba algorithm
solution = utils.solve(scrambled_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))