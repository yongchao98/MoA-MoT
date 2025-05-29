# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled state of the cube
scrambled_state = (
    'GGG', 'RYY', 'RYY',
    'YBBWRRGGWO', 'WRRGGWOOBB', 'WRRGGOYYYBB',
    'WWB', 'WWB', 'OOB'
)

# Convert the scrambled state into a format suitable for the solver
# The solver expects a single string with the colors in a specific order
# U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The given state needs to be mapped to this order
cube_string = (
    'GGGRYYRYY'  # U
    'WRRGGWOOBB'  # R
    'WRRGGOYYYBB'  # F
    'WWBOOB'  # D
    'YBBWRRGGWO'  # L
    'OOBWWB'  # B
)

# Solve the cube
solution = utils.solve(cube_string, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))