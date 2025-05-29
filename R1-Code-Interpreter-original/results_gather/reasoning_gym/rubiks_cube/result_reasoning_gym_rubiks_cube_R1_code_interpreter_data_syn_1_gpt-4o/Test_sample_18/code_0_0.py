# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube using the color notation
# The cube is represented as a string of 54 characters, each representing a color
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The given state is translated into this string format
scrambled_state = (
    "GRGGYGGYW"  # Up face
    "RRRWGBYBY"  # Right face
    "RRRWGWGOOY"  # Front face
    "BBRWWWGYBR"  # Down face
    "BBOOWOYBY"  # Left face
    "OWBOOWB"    # Back face
)

# Solve the cube using the solver
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to the required format (no combined steps)
formatted_solution = ' '.join(move.replace('2', ' ') for move in solution).replace('  ', ' ')

# Print the solution
print(formatted_solution)