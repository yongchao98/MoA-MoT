# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the cube configuration
# The cube is represented in a string format where each face is a single string
# The order is: U, R, F, D, L, B
# Each face is represented by 9 characters (3x3), and the whole cube is a single string of 54 characters
cube_configuration = (
    "YOR"  # U face
    "WYG"
    "RRG"
    "GRB"  # R face
    "WGO"
    "YYW"
    "GYR"  # F face
    "WRY"
    "BGW"
    "OOB"  # D face
    "YWB"
    "OOO"
    "WRY"  # L face
    "BGO"
    "OOW"
    "GYR"  # B face
    "BGB"
    "WRY"
)

# Solve the cube
solution = utils.solve(cube_configuration, 'Kociemba')

# Convert the solution to the required format
formatted_solution = ' '.join(solution).replace('2', ' ').replace("'", "' ")

# Print the solution
print(formatted_solution)