# Import the kociemba library
import kociemba

# Define the cube state in a string format
# The cube state is represented as a string of 54 characters, each representing a color on the cube
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# Colors: U=White, R=Red, F=Green, D=Yellow, L=Orange, B=Blue
# The given cube state is:
# U: G G G G Y Y W W B
# R: R R O G G W R B O
# F: Y Y Y W R O G G W
# D: W R O G G W R O Y
# L: W R B O O W R O Y
# B: R B B Y Y G B W W

# Convert the given cube state to the kociemba format
cube_state = (
    "GGGGYYWWB"  # U
    "RROGGWRBO"  # R
    "YYYWROGGW"  # F
    "WROGGWROY"  # D
    "WRBOOWROY"  # L
    "RBBYYGBWW"  # B
)

# Find the solution using the kociemba library
solution = kociemba.solve(cube_state)

# Print the solution
print(solution)