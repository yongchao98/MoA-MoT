# Importing a library to solve the Rubik's Cube
from rubik_solver import utils

# Define the scrambled cube state
# The cube state is represented as a string of 54 characters, each representing a color
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# Colors: W (White), R (Red), G (Green), Y (Yellow), O (Orange), B (Blue)
scrambled_state = (
    "BWG"  # Up face
    "RYO"
    "OOR"
    "OBG"  # Left face
    "WWW"
    "YBB"
    "OYY"  # Front face
    "RGY"
    "GOG"
    "WBR"  # Right face
    "RGY"
    "BOY"
    "WBW"  # Back face
    "RGY"
    "BWR"
    "BRR"  # Down face
)

# Solve the cube
solution = utils.solve(scrambled_state, 'Kociemba')

# Convert the solution to Singmaster notation without combining steps
singmaster_solution = ' '.join(move.replace('2', ' ') for move in solution)

# Print the solution
print(singmaster_solution)