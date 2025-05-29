# Initial board setup
board = [
    ['B', 'B', 'B', 'C', 'C', 'M'],
    ['.', '.', 'J', 'D', 'D', 'M'],
    ['A', 'A', 'J', 'K', '.', '.'],
    ['I', 'E', 'E', 'K', 'F', 'F'],
    ['I', 'G', 'G', 'K', 'L', 'x'],
    ['H', 'H', '.', '.', 'L', '.']
]

# Moves to solve the puzzle
moves = []

# Move M down
moves.append('M+1')

# Move C to the right
moves.append('C+1')

# Move B to the left
moves.append('B-1')

# Move A to the right to the exit
moves.append('A+3')

# Print the moves
print(' '.join(moves))