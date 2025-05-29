# Initial board setup
board = [
    ['B', 'B', 'H', '.', 'C', 'C'],
    ['.', 'G', 'H', '.', '.', 'K'],
    ['.', 'G', 'A', 'A', 'J', 'K'],
    ['D', 'D', '.', 'I', 'J', 'L'],
    ['.', '.', '.', 'I', 'x', 'L'],
    ['.', 'E', 'E', 'F', 'F', 'L']
]

# Moves to solve the puzzle
moves = []

# Move C to the right
moves.append('C+1')

# Move H down
moves.append('H+1')

# Move J down
moves.append('J+1')

# Move AA to the right to the exit
moves.append('A+3')

# Print the sequence of moves
print(' '.join(moves))