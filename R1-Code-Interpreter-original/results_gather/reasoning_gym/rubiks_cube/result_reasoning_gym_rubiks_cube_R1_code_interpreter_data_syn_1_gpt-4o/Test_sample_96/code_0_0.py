from rubik_solver import utils

# Define the current state of the cube using the given colors
# The state is represented as a string with each face's colors in order:
# U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The given cube state is:
# U: B Y B W Y B W G W
# R: Y R R B Y O G O O
# F: W O R G R R B G O
# D: Y B Y R W W Y G G
# L: G R R B G O G O O
# B: O Y O G Y G R B W

# Convert the given state into a string format for the solver
cube_state = (
    "BYBWYB"
    "WGWYRRBYOGOO"
    "WORGGRRBG"
    "YBYRWWYGG"
    "GRRBGOGOO"
    "OYO"
    "GYGRBW"
)

# Solve the cube using the solver library
solution = utils.solve(cube_state, 'Kociemba')

# Convert the solution to the required format (no combined steps)
formatted_solution = ' '.join(move.replace('2', ' ') for move in solution).replace('  ', ' ')

# Print the solution
print(formatted_solution)