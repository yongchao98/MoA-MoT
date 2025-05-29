# Importing the necessary library
from rubik_solver import utils

# Define the scrambled state of the cube
# The cube state is represented in a string format where each face is represented by a single letter
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The given cube state is:
# U: RRY YYY YYG
# R: OOB BBR GGO
# F: RRY RGG OOB
# D: WWW WWW WWB
# L: YOO BBB RRY
# B: GGO OOB RGG

# Convert the cube state into a string format
cube_state = "RRYYYYYYG" + "OOBBBRGGO" + "RRYRGGOOB" + "WWWWWWWWB" + "YOOBBBRRY" + "GGOOOBRGG"

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))