# Initial matrix
matrix = [
    [2, 0, 0],
    [0, 2, 0],
    [0, 16, 16]
]

# Initial position of the number to move
position = (2, 0)  # (row, column) in 0-indexed

# Move sequence
moves = 'dld'

# Simulate the moves
for move in moves:
    if move == 'd':
        # Move down
        if position[0] < 2:  # Check if not at the bottom
            position = (position[0] + 1, position[1])
    elif move == 'l':
        # Move left
        if position[1] > 0:  # Check if not at the leftmost
            position = (position[0], position[1] - 1)
    elif move == 'r':
        # Move right
        if position[1] < 2:  # Check if not at the rightmost
            position = (position[0], position[1] + 1)
    elif move == 'u':
        # Move up
        if position[0] > 0:  # Check if not at the top
            position = (position[0] - 1, position[1])

# Since the number is 0, the matrix remains unchanged
print(matrix)