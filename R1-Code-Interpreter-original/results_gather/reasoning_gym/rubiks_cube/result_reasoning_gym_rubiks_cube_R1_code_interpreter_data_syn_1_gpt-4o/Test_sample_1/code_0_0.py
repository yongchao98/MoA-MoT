from rubik_solver import utils

# Define the current state of the cube using the color notation
# The state is represented as a string of 54 characters, each representing a sticker
# The order is: UUUUUUUUURRRRRRRRRFFFFFFFFFDDDDDDDDDLLLLLLLLLBBBBBBBBB
# Using the given cube state:
# U: YYY YYY YYB
# R: RRR RRR GGR
# F: GGO GGO GGY
# D: WWW OWW OWW
# L: GGR GGR OOB
# B: OOO BBB RRB

cube_state = (
    "YYYYYYYYB"  # U
    "RRRRRRGGR"  # R
    "GGOGGOGGY"  # F
    "WWWOWWOWW"  # D
    "GGRGGROOB"  # L
    "OOOBBBRRB"  # B
)

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Print the solution in Singmaster notation
print(' '.join(solution))