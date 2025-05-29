# Define the initial board state
board = [
    "GBBIJL",
    "G.HIJL",
    "AAH.KL",
    "CCDDK.",
    "..x...",
    ".EEFF."
]

# Define the moves to solve the puzzle
moves = []

# Move H down by 1
moves.append("H+1")

# Move K down by 1
moves.append("K+1")

# Move L down by 1
moves.append("L+1")

# Move AA to the right by 4
moves.append("A+4")

# Print the moves in the required format
print("<<<" + " ".join(moves) + ">>>")