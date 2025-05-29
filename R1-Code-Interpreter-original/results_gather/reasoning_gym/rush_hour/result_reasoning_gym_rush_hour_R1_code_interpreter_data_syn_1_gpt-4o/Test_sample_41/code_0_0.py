# Define the initial board state
board = [
    "x.IBB.",
    "..IJCC",
    "AAIJK.",
    ".HDDKL",
    ".HEEEL",
    "FFGG.."
]

# Define the moves to solve the puzzle
moves = [
    "I-1",  # Move 'I' up by 1
    "J-1",  # Move 'J' up by 1
    "K-1",  # Move 'K' up by 1
    "L-1",  # Move 'L' up by 1
    "A+3"   # Move 'AA' to the right by 3
]

# Print the moves
print(" ".join(moves))