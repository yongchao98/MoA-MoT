# Define the initial board state
board = [
    "xH.BBB",
    ".H.xKL",
    "AA.JKL",
    "GCCJ.M",
    "G.IDDM",
    "EEIFFF"
]

# Define the moves to solve the puzzle
moves = [
    "J+1",  # Move J down
    "K+1",  # Move K down
    "L+1",  # Move L down
    "M+1",  # Move M down
    "A+2"   # Move AA to the right to the exit
]

# Print the moves
print(" ".join(moves))