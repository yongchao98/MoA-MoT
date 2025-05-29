# Initial board setup
board = [
    ['E', '.', '.', 'x', '.', '.'],
    ['E', 'B', 'B', 'G', '.', '.'],
    ['A', 'A', 'F', 'G', '.', '.'],
    ['.', '.', 'F', 'C', 'C', 'C'],
    ['.', '.', '.', '.', 'H', '.'],
    ['x', 'D', 'D', 'D', 'H', '.']
]

# Moves to solve the puzzle
moves = []

# Move G down by 1
moves.append('G+1')

# Move C right by 1
moves.append('C+1')

# Move F up by 1
moves.append('F-1')

# Move B up by 1
moves.append('B-1')

# Move E left by 1
moves.append('E-1')

# Move AA right by 4 to the exit
moves.append('A+4')

# Print the moves
print(' '.join(moves))