# Initial board setup
board = [
    "xBBHJK",
    "x..HJK",
    "FAAI.L",
    "FCCI.L",
    "..GDDL",
    "..GEEE"
]

# Moves to solve the puzzle
moves = []

# Move I down by 1
moves.append("I+1")

# Move L down by 1
moves.append("L+1")

# Move C up by 1
moves.append("C-1")

# Move AA to the right by 3
moves.append("A+3")

# Print the moves
print(" ".join(moves))