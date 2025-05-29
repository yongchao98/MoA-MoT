# Define the initial board state
board = [
    "HBBCCx",
    "H.IDDD",
    "AAIJ..",
    ".EEJFF",
    ".GGG.K",
    ".....K"
]

# Define the moves to solve the puzzle
moves = [
    "I+1",  # Move I down
    "J+1",  # Move J down
    "F+1",  # Move F right
    "A+3"   # Move A to the exit
]

# Print the moves
print(" ".join(moves))