# Initial board setup
board = [
    ['B', 'B', 'I', 'C', 'C', 'x'],
    ['.', '.', 'I', 'J', 'D', 'D'],
    ['A', 'A', 'I', 'J', '.', 'K'],
    ['H', 'E', 'E', 'F', 'F', 'K'],
    ['H', '.', 'G', 'G', '.', '.'],
    ['.', '.', '.', '.', '.', '.']
]

# Moves to solve the puzzle
moves = []

# Move K down
moves.append('K+1')

# Move D to the right
moves.append('D+1')

# Move A to the right to exit
moves.append('A+3')

# Print the moves
print(' '.join(moves))