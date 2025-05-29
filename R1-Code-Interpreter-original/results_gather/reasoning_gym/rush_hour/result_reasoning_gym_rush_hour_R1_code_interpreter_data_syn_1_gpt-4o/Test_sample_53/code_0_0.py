# Define the initial board state
board = [
    "..Ex..",
    "..EG.I",
    "AAFGHI",
    ".DF.H.",
    ".DBBH.",
    "..CCC."
]

# Define the moves to solve the puzzle
moves = [
    "F+1",  # Move F down
    "G+1",  # Move G down
    "H+1",  # Move H down
    "A+2"   # Move A to the right to the exit
]

# Print the moves
print(" ".join(moves))