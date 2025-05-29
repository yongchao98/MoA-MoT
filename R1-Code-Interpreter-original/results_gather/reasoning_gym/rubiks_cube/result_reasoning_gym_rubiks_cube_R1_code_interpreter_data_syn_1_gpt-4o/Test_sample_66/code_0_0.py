# Importing a library to solve the Rubik's Cube
from rubik_solver import utils

# Define the current state of the cube in a string format
# The cube state is represented in the order: U, R, F, D, L, B
# Each face is represented by 9 characters, starting from the top-left corner and going row by row
cube_state = (
    "BYOYYGGO"  # U face
    "YROYGGRRG"  # R face
    "OWGYGYBOW"  # F face
    "YGOBWWBRR"  # D face
    "RROYGYYOO"  # L face
    "GWRBBWBWB"  # B face
)

# Solve the cube using the library
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))