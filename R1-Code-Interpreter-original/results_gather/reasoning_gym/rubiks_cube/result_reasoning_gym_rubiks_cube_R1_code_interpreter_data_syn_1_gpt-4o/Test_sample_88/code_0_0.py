# Importing a Rubik's Cube solver library
from rubik_solver import utils

# Define the scrambled cube state
# The cube state is represented in a string format where each face is represented by a single letter
# The order of faces is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The given cube state is:
# U: WYYWYYWYY
# R: RRRRRRRRR
# F: BGGBGGBGG
# D: YWWYWWYWW
# L: OOOOOGGG
# B: BBBBBGGG

scrambled_state = 'WYYWYYWYYRRRRRRRRRBGGBGGBGGYWWYWWYWWOOOOOGGGBBBBGGG'

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))