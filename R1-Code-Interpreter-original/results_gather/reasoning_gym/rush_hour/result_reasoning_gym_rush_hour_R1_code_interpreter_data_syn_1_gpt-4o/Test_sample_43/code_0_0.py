# Initial board setup
board = [
    ['B', 'B', 'B', '.', 'x', 'K'],
    ['C', 'C', 'D', 'D', 'D', 'K'],
    ['H', 'A', 'A', 'J', '.', '.'],
    ['H', '.', '.', 'J', 'E', 'E'],
    ['.', '.', 'I', 'F', 'F', 'x'],
    ['G', 'G', 'I', '.', '.', '.']
]

# Moves to solve the puzzle
moves = []

# Move K down
moves.append('K+1')

# Move D left
moves.append('D-1')

# Move A to the right to exit
moves.append('A+3')

# Print the moves
print(' '.join(moves))