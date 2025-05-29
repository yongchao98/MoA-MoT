# Define the initial board state
board = [
    "..IBB.",
    ".GI.K.",
    ".GAAK.",
    "FGJCC.",
    "FHJ...",
    "FHDDEE"
]

# Define the moves
moves = []

# Move K down
moves.append("K+1")

# Move C to the right
moves.append("C+1")

# Move A to the right to exit
moves.append("A+2")

# Print the moves
print(" ".join(moves))