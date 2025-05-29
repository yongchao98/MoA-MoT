# Initial board setup
board = [
    ['x', 'B', 'B', 'B', '.', 'I'],
    ['.', 'C', 'C', '.', '.', 'I'],
    ['.', '.', 'G', 'A', 'A', 'I'],
    ['.', '.', 'G', 'H', 'D', 'D'],
    ['.', '.', 'G', 'H', 'E', 'E'],
    ['.', 'F', 'F', 'x', '.', '.']
]

# Moves to solve the puzzle
moves = []

# Move I down
moves.append('I+1')

# Move G up
moves.append('G-1')

# Move H up
moves.append('H-1')

# Move AA to the exit
moves.append('A+3')

# Print the moves
print(' '.join(moves))