# Importing the kociemba library to solve the Rubik's Cube
import kociemba

# The current state of the cube in a string format
# The cube state is represented in the order: U, R, F, D, L, B
# U: Up, R: Right, F: Front, D: Down, L: Left, B: Back
# Colors: W (White), Y (Yellow), R (Red), O (Orange), G (Green), B (Blue)
cube_state = (
    "YYGYYR"
    "WGB"
    "RROGOYOGOYBB"
    "RRYOGRWOOBBB"
    "BWWRWRWWGRGO"
    "GOBGWBWYY"
)

# Solving the cube using the kociemba library
solution = kociemba.solve(cube_state)

# Printing the solution in Singmaster notation
print(solution)