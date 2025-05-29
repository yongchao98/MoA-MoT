from rubik_solver import utils

# Define the current state of the cube
# The state is defined in the order: U, R, F, D, L, B
# Each face is represented as a string of 9 characters
cube_state = (
    "YYRGYYGWW",  # U
    "RRRRGGRGG",  # R
    "WBBWGBOOB",  # F
    "WWYWWYRB",   # D
    "GGGOWOOGO",  # L
    "OBBWROOY"    # B
)

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))