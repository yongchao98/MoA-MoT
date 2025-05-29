# Initial board setup
board = [
    ['.', '.', '.', '.', 'B', 'B'],
    ['.', '.', '.', 'x', '.', 'F'],
    ['A', 'A', '.', 'E', '.', 'F'],
    ['.', '.', '.', 'E', '.', '.'],
    ['C', 'C', 'D', 'D', '.', '.'],
    ['.', '.', '.', '.', '.', '.']
]

# Moves to solve the puzzle
moves = []

# Move E down by 1
moves.append('E+1')

# Move F to the right by 1
moves.append('F+1')

# Move AA to the right by 2 to reach the exit
moves.append('A+2')

# Print the moves
print(' '.join(moves))