# Initial board setup
board = [
    "GBB.KL",
    "GHI.KL",
    "GHIAAM",
    "CCCJ.M",
    "..xJDD",
    "EEFF.."
]

# Moves to solve the puzzle
moves = [
    "M+2",  # Move M down by 2 spaces
    "J+1",  # Move J down by 1 space
    "C-3",  # Move C left by 3 spaces
    "A+1"   # Move A right by 1 space to the exit
]

# Print the moves
print(" ".join(moves))