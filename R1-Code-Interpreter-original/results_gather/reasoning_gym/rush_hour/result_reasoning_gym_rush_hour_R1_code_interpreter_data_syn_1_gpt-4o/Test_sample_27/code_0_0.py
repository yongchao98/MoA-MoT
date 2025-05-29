# Initial board setup
board = [
    ['F', 'B', 'B', '.', 'J', '.'],
    ['F', '.', 'I', '.', 'J', 'K'],
    ['A', 'A', 'I', '.', '.', 'K'],
    ['C', 'C', 'I', '.', '.', 'L'],
    ['G', 'H', 'D', 'D', 'D', 'L'],
    ['G', 'H', 'E', 'E', 'x', '.']
]

# Moves to solve the puzzle
moves = []

# Move 'I' down
moves.append('I+1')

# Move 'J' to the right
moves.append('J+1')

# Move 'K' to the right
moves.append('K+1')

# Move 'AA' to the right
moves.append('A+4')

# Print the moves
print(' '.join(moves))