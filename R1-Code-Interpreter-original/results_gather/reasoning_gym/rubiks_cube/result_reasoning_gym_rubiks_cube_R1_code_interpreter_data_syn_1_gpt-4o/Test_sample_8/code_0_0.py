# Import the kociemba library
import kociemba

# Define the current state of the cube
# The cube state is represented as a string of 54 characters
# Each character represents a sticker on the cube
# The order is: U (up), R (right), F (front), D (down), L (left), B (back)
cube_state = (
    "YRRRYWGGR"
    "BWRWGBWOYBYO"
    "ORYOGBWOWGBG"
    "GGO"
    "BWR"
    "WOG"
)

# Use the kociemba library to find the solution
solution = kociemba.solve(cube_state)

# Print the solution in Singmaster notation
print(solution)