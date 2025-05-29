# Initial matrix
matrix = [
    [8, 32],
    [0, 2]
]

# Move sequence
moves = 'uu'

# Initial position of the number to move
position = (1, 1)  # (row, column) in zero-indexed

# Perform the moves
for move in moves:
    if move == 'u' and position[0] > 0:
        # Move up
        matrix[position[0] - 1][position[1]] = matrix[position[0]][position[1]]
        matrix[position[0]][position[1]] = 0
        position = (position[0] - 1, position[1])

# Print the final matrix
print(matrix)