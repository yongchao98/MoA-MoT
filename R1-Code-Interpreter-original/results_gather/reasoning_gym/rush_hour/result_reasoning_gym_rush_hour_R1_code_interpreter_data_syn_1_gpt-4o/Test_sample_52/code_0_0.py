# Define the initial board state
board = [
    "BBH...",
    "x.HCCC",
    ".AAI.K",
    ".G.IJK",
    ".GDDJK",
    ".EEEFF"
]

# Define the moves
moves = []

# Move K down by 1
moves.append("K+1")

# Move J down by 1
moves.append("J+1")

# Move A right by 2 to exit
moves.append("A+2")

# Print the moves
print(" ".join(moves))