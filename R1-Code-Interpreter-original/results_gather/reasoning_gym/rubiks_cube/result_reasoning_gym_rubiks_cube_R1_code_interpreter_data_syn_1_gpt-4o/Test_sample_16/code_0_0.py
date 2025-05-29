from rubik_solver import utils

# Define the scrambled cube state in a string format
# The cube state is represented in the order: U, R, F, D, L, B
# Each face is represented by 9 characters, starting from the top-left corner and going row by row
scrambled_state = (
    "BYYGYYRYG"  # U face
    "GOWYRRGGY"  # R face
    "YRRGGYWOW"  # F face
    "WOOWWRRB"   # D face
    "YRRGOOWOB"  # L face
    "WOBBWRRBB"  # B face
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))