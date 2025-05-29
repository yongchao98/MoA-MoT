# Initial board setup
board = [
    ['.', '.', '.', '.', 'H', '.'],
    ['.', 'x', 'B', 'B', 'H', '.'],
    ['.', 'A', 'A', '.', 'I', '.'],
    ['C', 'C', 'F', 'G', 'I', '.'],
    ['E', '.', 'F', 'G', 'D', 'D'],
    ['E', '.', 'x', '.', '.', '.']
]

# Moves to solve the puzzle
moves = []

# Move H to the right
moves.append('H+1')

# Move B up
moves.append('B-1')

# Move A to the right to the exit
moves.append('A+3')

# Print the moves
print(' '.join(moves))