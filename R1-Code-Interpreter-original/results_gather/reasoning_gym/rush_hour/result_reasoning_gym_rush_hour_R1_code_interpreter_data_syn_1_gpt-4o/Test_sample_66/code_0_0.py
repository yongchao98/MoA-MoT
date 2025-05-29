# Define the initial board state
board = [
    ['x', '.', 'G', 'B', 'B', 'J'],
    ['C', 'C', 'G', 'H', 'I', 'J'],
    ['F', 'A', 'A', 'H', 'I', 'K'],
    ['F', 'D', 'D', '.', 'I', 'K'],
    ['E', 'E', '.', 'x', '.', '.'],
    ['.', '.', '.', '.', '.', '.']
]

# Define the moves to solve the puzzle
moves = []

# Move J down
moves.append('J+1')

# Move G left
moves.append('G-1')

# Move AA to the right
moves.append('A+3')

# Print the moves
print(' '.join(moves))