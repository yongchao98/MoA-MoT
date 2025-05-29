# Initial matrix
matrix = [
    [0, 16, 32],
    [0, 0, 32],
    [16, 4, 2]
]

# Initial position of the number
row, col = 2, 1  # 0-indexed position (3, 2)

# Move sequence
moves = 'dlr'

# Execute the moves
for move in moves:
    if move == 'd':
        # No movement as it's already at the bottom
        continue
    elif move == 'l':
        if col > 0:  # Check if we can move left
            # Move the number and check for combination
            if matrix[row][col - 1] == matrix[row][col]:
                matrix[row][col - 1] *= 2
                matrix[row][col] = 0
            else:
                matrix[row][col - 1], matrix[row][col] = matrix[row][col], matrix[row][col - 1]
            col -= 1
    elif move == 'r':
        if col < 2:  # Check if we can move right
            # Move the number and check for combination
            if matrix[row][col + 1] == matrix[row][col]:
                matrix[row][col + 1] *= 2
                matrix[row][col] = 0
            else:
                matrix[row][col + 1], matrix[row][col] = matrix[row][col], matrix[row][col + 1]
            col += 1

# Print the final matrix
print(matrix)