# Import the kociemba library
import kociemba

# Define the scrambled state of the cube
# The state is represented as a string with each face's colors in order:
# U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The colors are represented by their initials: W (White), R (Red), B (Blue), O (Orange), G (Green), Y (Yellow)

scrambled_state = (
    "OBBYYYYYY"
    "WBBRRRGGO"
    "RRRGGYOOG"
    "WWGWWGWWR"
    "RRRGGYOOG"
    "BBBRRYOOO"
)

# Solve the cube using the kociemba library
solution = kociemba.solve(scrambled_state)

# Print the solution
print(solution)