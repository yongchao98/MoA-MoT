# Import the kociemba library
import kociemba

# Define the scrambled cube state
# The cube state is represented as a string of 54 characters, each representing a color on the cube.
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# Colors: U=Y, R=R, F=G, D=W, L=O, B=B

scrambled_state = (
    "YYYYYWB"  # Up face
    "RRRGGWRROOBBB"  # Right, Front, Left faces
    "OOWGWRROOBBY"  # Down face
    "GYGGWWGWW"  # Back face
)

# Find the solution using the kociemba library
solution = kociemba.solve(scrambled_state)

# Print the solution
print(solution)