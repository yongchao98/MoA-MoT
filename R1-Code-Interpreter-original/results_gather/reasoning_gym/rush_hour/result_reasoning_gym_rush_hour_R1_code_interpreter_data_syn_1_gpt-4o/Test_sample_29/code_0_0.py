# Initial board setup
board = [
    ['B', 'B', 'B', 'C', 'C', '.'],
    ['.', '.', 'J', 'K', 'D', 'D'],
    ['A', 'A', 'J', 'K', 'L', 'M'],
    ['I', '.', 'E', 'E', 'L', 'M'],
    ['I', 'F', 'F', '.', '.', 'N'],
    ['G', 'G', 'H', 'H', 'H', 'N']
]

# Moves to solve the puzzle
moves = []

# Move D down
moves.append('D+1')

# Move L down
moves.append('L+1')

# Move M down
moves.append('M+1')

# Move A to the exit
moves.append('A+4')

# Print the sequence of moves
print(' '.join(moves))