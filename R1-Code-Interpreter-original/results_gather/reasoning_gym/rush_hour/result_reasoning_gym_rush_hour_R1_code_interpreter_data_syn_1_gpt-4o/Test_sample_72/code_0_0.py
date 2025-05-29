# Initial board setup
board = [
    ['.', '.', 'x', 'B', 'B', 'x'],
    ['.', '.', '.', 'E', '.', '.'],
    ['A', 'A', '.', 'E', '.', '.'],
    ['.', 'C', 'C', 'C', '.', '.'],
    ['.', '.', 'D', 'D', '.', '.'],
    ['.', '.', '.', '.', '.', '.']
]

# Moves to solve the puzzle
moves = []

# Move E down by 1
moves.append('E+1')

# Move AA right by 2
moves.append('A+2')

# Print the moves
print(' '.join(moves))