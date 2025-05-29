# Define the initial board state
board = [
    "..BB..",
    "xCC.x.",
    "AAG...",
    "F.GDD.",
    "FEEH..",
    "F..H.."
]

# Define the moves to solve the puzzle
moves = []

# Move D left by 2 spaces
moves.append("D-2")

# Move G down by 1 space
moves.append("G+1")

# Move A (red car) right by 3 spaces to the exit
moves.append("A+3")

# Print the moves
print(" ".join(moves))