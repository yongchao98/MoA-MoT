# Initial board setup
board = [
    "...HI.",
    "BBBHI.",
    "AA.HJ.",
    "F.CCJ.",
    "F.GDD.",
    "EEG..."
]

# Moves to solve the puzzle
moves = []

# Step 1: Move 'J' down by 1
moves.append('J+1')

# Step 2: Move 'C' down by 1
moves.append('C+1')

# Step 3: Move 'A' to the right by 2 to reach the exit
moves.append('A+2')

# Print the sequence of moves
print(' '.join(moves))