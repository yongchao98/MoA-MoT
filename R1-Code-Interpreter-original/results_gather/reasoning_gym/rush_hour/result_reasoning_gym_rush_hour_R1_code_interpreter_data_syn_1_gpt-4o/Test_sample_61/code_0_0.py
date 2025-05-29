# Initial board setup
board = [
    ['B', 'B', 'B', 'C', 'C', 'M'],
    ['D', 'D', 'E', 'E', 'E', 'M'],
    ['A', 'A', 'J', 'K', '.', '.'],
    ['I', '.', 'J', 'K', 'F', 'F'],
    ['I', '.', 'G', 'G', 'L', 'x'],
    ['H', 'H', '.', '.', 'L', '.']
]

# Moves to solve the puzzle
moves = []

# Move M down by 1
moves.append('M+1')

# Move E down by 1
moves.append('E+1')

# Move M down by 1 again
moves.append('M+1')

# Move AA right by 3 to exit
moves.append('A+3')

# Print the moves
print(' '.join(moves))