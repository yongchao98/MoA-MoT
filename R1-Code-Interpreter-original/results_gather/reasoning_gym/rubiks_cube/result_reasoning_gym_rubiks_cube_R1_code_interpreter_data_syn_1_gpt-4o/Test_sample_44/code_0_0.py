# Import the kociemba library
import kociemba

# Define the scrambled state of the cube
# The state is represented as a string of 54 characters, each representing a color on the cube
# The order is: U (Up), R (Right), F (Front), D (Down), L (Left), B (Back)
# The colors are represented by their initials: W (White), Y (Yellow), R (Red), O (Orange), G (Green), B (Blue)

scrambled_state = (
    "OGRWYYGWW"  # U face
    "RYOOGOYBR"  # R face
    "RGWGBGOOB"  # F face
    "RYRBWYOYW"  # D face
    "BROWRBOWY"  # L face
    "GOYBRWGBY"  # B face
)

# Find the solution using the kociemba library
solution = kociemba.solve(scrambled_state)

# Print the solution
print(solution)