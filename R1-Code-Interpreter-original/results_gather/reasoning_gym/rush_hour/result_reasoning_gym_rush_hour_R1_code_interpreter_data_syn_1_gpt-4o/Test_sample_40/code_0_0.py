# Define the initial board state
board = [
    ['.', 'E', 'B', 'B', 'H', '.'],
    ['.', 'E', '.', 'G', 'H', '.'],
    ['A', 'A', 'G', '.', '.', '.'],
    ['.', '.', 'F', 'C', 'C', 'C'],
    ['.', '.', 'F', '.', '.', '.'],
    ['.', 'D', 'D', 'D', '.', 'x']
]

# Define the moves to solve the puzzle
moves = [
    'G+1',  # Move G down
    'C+1',  # Move C to the right
    'A+3'   # Move A to the right to exit
]

# Print the moves
print(' '.join(moves))