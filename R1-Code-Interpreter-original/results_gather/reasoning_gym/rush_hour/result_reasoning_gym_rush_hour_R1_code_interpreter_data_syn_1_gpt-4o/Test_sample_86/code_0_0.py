# Initial board setup
board = [
    ['.', 'B', 'B', '.', 'x', 'M'],
    ['.', '.', 'J', 'C', 'C', 'M'],
    ['A', 'A', 'J', '.', 'L', 'M'],
    ['D', 'D', '.', 'K', 'L', '.'],
    ['I', 'E', 'E', 'K', 'F', 'F'],
    ['I', 'G', 'G', 'G', 'H', 'H']
]

# Moves to solve the puzzle
moves = []

# Move J down
moves.append('J+1')

# Move M to the right
moves.append('M+1')

# Move AA to the right to the exit
moves.append('A+4')

# Print the sequence of moves
print(' '.join(moves))