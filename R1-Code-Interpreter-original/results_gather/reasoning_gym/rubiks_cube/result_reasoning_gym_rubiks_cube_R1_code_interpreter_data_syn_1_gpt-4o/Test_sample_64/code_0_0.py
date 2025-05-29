# Importing a library to solve the Rubik's Cube
from kociemba import solve

# Define the cube state in a string format
# The cube state is represented as a string of 54 characters, each representing a color
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# Colors: W (White), R (Red), B (Blue), O (Orange), G (Green), Y (Yellow)

cube_state = (
    "GYBGYYGOG"
    "YRRWGRGOW"
    "RRRWGWBOW"
    "BGGOYOBOR"
    "WGW"
    "WWBOBB"
)

# Solve the cube using the kociemba library
solution = solve(cube_state)

# Print the solution in Singmaster notation
print(solution)