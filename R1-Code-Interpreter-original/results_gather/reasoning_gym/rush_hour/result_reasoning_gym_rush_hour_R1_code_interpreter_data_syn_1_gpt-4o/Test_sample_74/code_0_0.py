# Initial board setup
board = [
    ['B', 'B', 'B', 'C', 'C', '.'],
    ['.', '.', 'H', 'I', 'D', 'D'],
    ['A', 'A', 'H', 'I', 'J', 'K'],
    ['G', 'E', 'E', 'E', 'J', 'K'],
    ['G', '.', 'F', 'F', '.', 'x'],
    ['.', '.', '.', '.', '.', '.']
]

# Moves to solve the puzzle
moves = []

# Move 'H' down by 1
moves.append('H+1')

# Move 'I' down by 1
moves.append('I+1')

# Move 'J' down by 1
moves.append('J+1')

# Move 'AA' right by 4
moves.append('A+4')

# Print the moves
print(' '.join(moves))