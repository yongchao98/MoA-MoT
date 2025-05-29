# Initial matrix
matrix = [
    [8, 0],
    [4, 0]
]

# Move sequence
sequence = 'ur'

# Initial position of the number to move
position = (0, 0)  # (row, column) in 0-indexed

# Execute the moves
for move in sequence:
    if move == 'u':
        # Move up: decrease row index if possible
        if position[0] > 0:
            position = (position[0] - 1, position[1])
    elif move == 'r':
        # Move right: increase column index if possible
        if position[1] < 1:
            position = (position[0], position[1] + 1)

# Update the matrix with the new position
new_matrix = [[0, 0], [0, 0]]
new_matrix[position[0]][position[1]] = matrix[0][0]  # Move the number 8 to the new position
new_matrix[1][0] = matrix[1][0]  # Keep the number 4 in its original position

print(new_matrix)