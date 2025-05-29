# Define the initial board state
board = [
    "GBBIJK",
    "G..IJK",
    "AAHI..",
    "..HCCC",
    "..xDD.",
    "EEEFF."
]

# Define the moves
moves = []

# Move 'H' down by 1
moves.append("H+1")

# Move 'C' left by 1
moves.append("C-1")

# Move 'A' right by 3 to exit
moves.append("A+3")

# Print the moves in the required format
print("<<<" + " ".join(moves) + ">>>")