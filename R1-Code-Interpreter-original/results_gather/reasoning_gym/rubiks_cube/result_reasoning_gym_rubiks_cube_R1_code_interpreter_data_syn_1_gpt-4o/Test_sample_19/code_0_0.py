# Importing a library to solve the Rubik's Cube
from kociemba import solve

# The cube's current state in a string format
# The string is constructed by reading the cube's faces in the order: U, R, F, D, L, B
# Each face is read in row-major order
cube_state = (
    "BRBYWYWGOWBYRWROBGYOBORRYWWBGOGOWRYYRGGWOGYOBY"
)

# Solving the cube
solution = solve(cube_state)

# Printing the solution in Singmaster notation
print(solution)