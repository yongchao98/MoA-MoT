# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the cube configuration in a format that the solver can understand
# The cube is represented as a string of 54 characters, each representing a color
# The order is: U (up), R (right), F (front), D (down), L (left), B (back)
# The colors are represented by their initials: W, R, B, G, Y, O

cube_configuration = (
    "WGG"  # U face
    "YYO"
    "YYO"
    "OOOBBGWWRWWB"  # L, F, R faces
    "RRYRGBOOOGBB"
    "RRYRGOGGGRBB"
    "BRY"  # D face
    "WWY"
    "WWY"
)

# Solve the cube using the solver
solution = utils.solve(cube_configuration, 'Kociemba')

# Convert the solution to the required format (Singmaster notation without combining steps)
singmaster_solution = ' '.join(move.replace('2', ' 2') for move in solution)

# Print the solution
print(singmaster_solution)