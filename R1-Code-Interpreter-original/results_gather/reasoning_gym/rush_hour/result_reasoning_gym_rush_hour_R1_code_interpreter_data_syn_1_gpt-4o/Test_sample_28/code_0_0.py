# Define the initial board state
board = [
    "..xBBJ",
    "CCCDDJ",
    "..AA.K",
    "G.HEEK",
    "G.HIFF",
    "...I.."
]

# Define the moves
moves = []

# Move 'K' down by 1
moves.append('K+1')

# Move 'AA' right by 3 to exit
moves.append('A+3')

# Print the moves
print(' '.join(moves))