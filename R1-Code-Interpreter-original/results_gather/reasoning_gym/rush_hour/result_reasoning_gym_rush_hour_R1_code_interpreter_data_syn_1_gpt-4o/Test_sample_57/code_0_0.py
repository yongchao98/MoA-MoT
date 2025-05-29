# Define the initial board state
board = [
    ".xEBBB",
    "..EFG.",
    "AA.FGH",
    "D.CC.H",
    "D....x",
    "D....."
]

# Define the moves to solve the puzzle
moves = []

# Move F down by 1
moves.append("F+1")

# Move G down by 1
moves.append("G+1")

# Move H down by 1
moves.append("H+1")

# Move A to the right by 3 to exit
moves.append("A+3")

# Print the moves in the required format
print("<<<" + " ".join(moves) + ">>>")