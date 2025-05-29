# Define the initial board state
board = [
    ['.', 'F', 'B', 'B', 'J', 'K'],
    ['.', 'F', '.', 'I', 'J', 'K'],
    ['x', 'A', 'A', 'I', '.', 'L'],
    ['C', 'C', 'H', 'I', '.', 'L'],
    ['E', 'G', 'H', '.', '.', '.'],
    ['E', 'G', 'D', 'D', 'D', '.']
]

# Define the moves
moves = []

# Move 'I' down by 1
moves.append('I+1')

# Move 'AA' right by 3
moves.append('A+3')

# Print the moves
print(' '.join(moves))