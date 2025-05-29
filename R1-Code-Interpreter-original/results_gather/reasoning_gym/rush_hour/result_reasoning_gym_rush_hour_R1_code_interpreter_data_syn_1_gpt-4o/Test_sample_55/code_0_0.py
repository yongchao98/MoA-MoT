# Initial board setup
board = [
    ['E', 'B', 'B', '.', 'C', 'C'],
    ['E', '.', 'F', '.', '.', 'H'],
    ['A', 'A', 'F', '.', 'G', 'H'],
    ['.', '.', '.', '.', 'G', 'I'],
    ['.', 'x', 'D', 'D', 'G', 'I'],
    ['.', '.', '.', '.', '.', 'I']
]

# Moves to solve the puzzle
moves = []

# Move F down
moves.append('F+1')

# Move G down
moves.append('G+1')

# Move H down
moves.append('H+1')

# Move I down
moves.append('I+1')

# Move A to the right
moves.append('A+3')

# Print the sequence of moves
print(' '.join(moves))