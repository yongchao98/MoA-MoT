# Define the initial board state
board = [
    ".....x",
    "BB.xJ.",
    "FAAIJ.",
    "F.HICC",
    "G.HDDK",
    "GEEE.K"
]

# Define the moves to solve the puzzle
moves = [
    "J+1",  # Move J down
    "I+1",  # Move I down
    "C+1",  # Move C down
    "A+2"   # Move AA to the right to the exit
]

# Print the moves in the required format
print("<<<" + " ".join(moves) + ">>>")