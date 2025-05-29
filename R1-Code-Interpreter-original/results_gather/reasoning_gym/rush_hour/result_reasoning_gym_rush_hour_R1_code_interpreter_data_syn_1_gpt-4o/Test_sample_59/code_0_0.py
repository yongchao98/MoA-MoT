# Initial board setup
board = [
    ['B', 'B', '.', '.', '.', 'J'],
    ['.', 'G', '.', '.', 'I', 'J'],
    ['.', 'G', 'A', 'A', 'I', 'J'],
    ['.', '.', '.', 'H', 'C', 'C'],
    ['D', 'D', '.', 'H', '.', '.'],
    ['E', 'E', 'F', 'F', 'x', '.']
]

# Moves to solve the puzzle
moves = []

# Step 1: Move 'J' down
moves.append('J+1')

# Step 2: Move 'I' to the right
moves.append('I+1')

# Step 3: Move 'AA' to the right to the exit
moves.append('A+3')

# Print the moves
print(' '.join(moves))