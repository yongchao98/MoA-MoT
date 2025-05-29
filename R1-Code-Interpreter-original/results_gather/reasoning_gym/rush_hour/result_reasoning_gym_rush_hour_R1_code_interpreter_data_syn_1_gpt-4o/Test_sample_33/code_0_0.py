# Initial board setup
board = [
    ['.', '.', 'x', 'B', 'B', 'K'],
    ['G', 'C', 'C', 'J', '.', 'K'],
    ['G', 'A', 'A', 'J', '.', '.'],
    ['D', 'D', 'I', 'J', '.', '.'],
    ['.', 'H', 'I', 'E', 'E', '.'],
    ['.', 'H', 'F', 'F', 'x', '.']
]

# Moves to solve the puzzle
moves = []

# Move J down
moves.append('J+1')

# Move K to the right
moves.append('K+1')

# Move AA to the right to the exit
moves.append('A+3')

# Print the moves
print(' '.join(moves))