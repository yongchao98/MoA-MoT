# Initial board setup
board = [
    ['G', 'B', 'B', '.', 'L', 'M'],
    ['G', 'H', 'I', '.', 'L', 'M'],
    ['G', 'H', 'I', 'A', 'A', 'N'],
    ['C', 'C', 'J', 'K', '.', 'N'],
    ['.', '.', 'J', 'K', 'D', 'D'],
    ['.', 'E', 'E', 'F', 'F', '.']
]

# Moves to solve the puzzle
moves = []

# Move N down
moves.append('N+1')

# Move M to the right
moves.append('M+1')

# Move A to the right to the exit
moves.append('A+1')

# Print the moves
print(' '.join(moves))