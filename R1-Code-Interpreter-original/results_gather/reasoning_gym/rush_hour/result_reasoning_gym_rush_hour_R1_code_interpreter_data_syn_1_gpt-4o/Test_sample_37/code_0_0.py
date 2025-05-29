# Define the initial board state
board = [
    ['B', 'B', 'C', 'C', 'C', '.'],
    ['.', '.', '.', 'G', 'H', '.'],
    ['A', 'A', '.', 'G', 'H', '.'],
    ['.', '.', 'D', 'D', 'H', '.'],
    ['.', 'E', 'E', 'E', '.', '.'],
    ['.', '.', 'F', 'F', 'F', '.']
]

# Define the moves to solve the puzzle
moves = [
    'B-1',  # Move 'BB' up
    'C+1',  # Move 'CCC' right
    'D+1',  # Move 'DD' down
    'E-1',  # Move 'EEE' up
    'F+1',  # Move 'FFF' down
    'G+1',  # Move 'GH' right
    'A+4'   # Move 'AA' to the exit
]

# Print the moves
print(' '.join(moves))