# Initial board setup
board = [
    "xBBHJK",
    "FCCHJK",
    "FAAI.L",
    "F..I.L",
    "..GDDL",
    "EEG..."
]

# Function to find the moves
def find_moves():
    moves = []
    # Move 'I' down by 1
    moves.append('I+1')
    # Move 'L' down by 1
    moves.append('L+1')
    # Move 'AA' to the right by 3
    moves.append('A+3')
    return moves

# Get the moves
moves = find_moves()

# Print the moves in the required format
print(f"<<<{' '.join(moves)}>>>")