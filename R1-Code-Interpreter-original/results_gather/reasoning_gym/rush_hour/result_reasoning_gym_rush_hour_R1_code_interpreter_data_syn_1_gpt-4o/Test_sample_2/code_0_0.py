# Initial board setup
board = [
    ['B', 'B', 'H', '.', 'C', 'C'],
    ['G', '.', 'H', '.', 'J', 'K'],
    ['G', 'A', 'A', 'J', 'K'],
    ['D', 'D', '.', 'I', 'x', 'L'],
    ['E', 'E', '.', 'I', '.', 'L'],
    ['F', 'F', 'F', '.', '.', 'x']
]

# Moves to solve the puzzle
moves = []

# Move J down
moves.append('J+1')

# Move K down
moves.append('K+1')

# Move A (red car) to the right
moves.append('A+3')

# Print the moves
print(' '.join(moves))