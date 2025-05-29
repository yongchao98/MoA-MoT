# Import the kociemba library
import kociemba

# Define the current state of the cube in the kociemba format
# The format is a string of 54 characters representing the cube's faces in the order:
# UUUUUUUUURRRRRRRRRFFFFFFFFFDDDDDDDDDLLLLLLLLLBBBBBBBBB
# Each face is represented by its initial: U, R, F, D, L, B
# The given cube state needs to be converted to this format

# Current state of the cube
cube_state = (
    "ORYYYYBBB"  # U face
    "BBWRRRYGG"  # R face
    "OWYOWGOBG"  # F face
    "BWGBWGRGW"  # D face
    "WRRYGWOOG"  # L face
    "OOWYBWRRY"  # B face
)

# Find the solution using the kociemba library
solution = kociemba.solve(cube_state)

# Print the solution
print(solution)