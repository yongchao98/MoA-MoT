# Initial board setup
board = [
    "..HBBL",
    "..HI.L",
    "AAHIKL",
    "GCCCK.",
    "GDDJ..",
    "EE.JFF"
]

# Moves to solve the puzzle
moves = []

# Move L down
moves.append("L+1")

# Move K down
moves.append("K+1")

# Move C left
moves.append("C-1")

# Move A to the exit
moves.append("A+4")

# Print the moves
print(" ".join(moves))