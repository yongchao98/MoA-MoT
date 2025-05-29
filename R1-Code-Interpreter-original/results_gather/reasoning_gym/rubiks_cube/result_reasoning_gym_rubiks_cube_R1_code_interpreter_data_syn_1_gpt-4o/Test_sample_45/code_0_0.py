from rubik_solver import utils

# Define the scrambled state of the cube
# The state is defined in the order: U, R, F, D, L, B
# Each face is represented by a string of 9 characters
# U: Upper face, R: Right face, F: Front face, D: Down face, L: Left face, B: Back face
scrambled_state = (
    "RYRGYGWWW",  # U
    "OOGGRBWRR",  # R
    "GYYRRBWRB",  # F
    "BBOWWW",     # D
    "GOGRRBOOY",  # L
    "GBYBOOY"     # B
)

# Solve the cube
solution = utils.solve(scrambled_state)

# Print the solution in Singmaster notation
print(' '.join(solution))