# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the current state of the cube in a format that the solver can understand
# The format is a string of 54 characters representing the cube's faces in the order:
# UUUUUUUUURRRRRRRRRFFFFFFFFFDDDDDDDDDLLLLLLLLLBBBBBBBBB
# Each face is represented by 9 characters, and the order is Up, Right, Front, Down, Left, Back

# Current state of the cube
cube_state = (
    "RWGBYRBRBWBYGRGYGYBOYWRBWBYWWWGRBWOGOWGOWGO"
)

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))