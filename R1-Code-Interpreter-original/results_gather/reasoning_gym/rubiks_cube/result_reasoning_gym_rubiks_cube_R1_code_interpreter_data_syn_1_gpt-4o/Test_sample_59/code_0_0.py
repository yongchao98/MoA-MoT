# Import the kociemba library
import kociemba

# Define the cube state in the kociemba format
# The cube state is represented as a string of 54 characters, each representing a color
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The colors are represented by their initials: W (White), R (Red), B (Blue), O (Orange), G (Green), Y (Yellow)

# Given cube state:
# U: O Y R O Y R G Y Y
# R: G G O G G Y O W W
# F: Y R Y R G G O O W
# D: B W W W W W R R R
# L: B B R W G G O G G
# B: Y R O Y G G O O W

# Convert the given state to the kociemba format
cube_state = (
    "OYR"  # U face
    "OYR"
    "GYY"
    "GGG"  # R face
    "GGY"
    "OWW"
    "YRY"  # F face
    "RGG"
    "OOW"
    "BWW"  # D face
    "WWW"
    "RRR"
    "BBR"  # L face
    "WGG"
    "OGG"
    "YRO"  # B face
    "YGG"
    "OOW"
)

# Find the solution using the kociemba library
solution = kociemba.solve(cube_state)

# Print the solution
print(solution)