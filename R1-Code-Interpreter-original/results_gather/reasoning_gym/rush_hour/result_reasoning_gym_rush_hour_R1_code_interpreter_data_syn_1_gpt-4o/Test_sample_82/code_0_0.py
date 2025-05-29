# Initial board setup
board = [
    ['F', '.', '.', 'H', 'B', 'B'],
    ['F', 'C', 'C', 'H', '.', 'I'],
    ['A', 'A', 'G', '.', '.', 'I'],
    ['.', '.', 'G', 'D', 'D', 'I'],
    ['.', '.', 'x', 'x', '.', '.'],
    ['.', '.', '.', 'E', 'E', 'E']
]

# Moves to solve the puzzle
moves = []

# Move H down by 1
moves.append('H+1')

# Move C left by 2
moves.append('C-2')

# Move A right by 3 to the exit
moves.append('A+3')

# Print the moves
print(' '.join(moves))