# Initial board setup
board = [
    ['B', 'B', 'H', '.', 'J', '.'],
    ['.', '.', 'H', 'x', 'J', 'K'],
    ['A', 'A', 'H', '.', '.', 'K'],
    ['G', 'C', 'C', 'C', '.', 'K'],
    ['G', 'D', 'D', 'I', '.', '.'],
    ['E', 'E', '.', 'I', 'F', 'F']
]

# Moves to solve the puzzle
moves = []

# Move J down
moves.append('J+1')

# Move H up
moves.append('H-1')

# Move AA to the exit
moves.append('A+4')

# Print the moves
print(' '.join(moves))