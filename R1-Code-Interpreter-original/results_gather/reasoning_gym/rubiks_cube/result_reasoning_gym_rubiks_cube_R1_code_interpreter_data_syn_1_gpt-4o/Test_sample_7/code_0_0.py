from rubik_solver import utils

# Define the cube state in a format that the solver can understand
# The cube state is represented as a string of 54 characters, each representing a color
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The given cube state is:
# U: R Y R B Y B O O O
# R: G R B Y Y Y G O B
# F: W B W G R Y G G O
# D: R R R G W G O W O
# L: G R B Y Y Y G O B
# B: W B W W B W W R W

# Convert the given cube state to a single string
cube_state = (
    "RYRBYBOOO"  # U
    "GRBYYYGOB"  # R
    "WBWGRYGGO"  # F
    "RRRGWGOWO"  # D
    "GRBYYYGOB"  # L
    "WBWWBWWRW"  # B
)

# Solve the cube
solution = utils.solve(cube_state, 'Kociemba')

# Convert the solution to the required format (no combined moves)
formatted_solution = ' '.join(move.replace('2', ' ').replace('\'', ' \'') for move in solution).replace('  ', ' ')

print(formatted_solution)