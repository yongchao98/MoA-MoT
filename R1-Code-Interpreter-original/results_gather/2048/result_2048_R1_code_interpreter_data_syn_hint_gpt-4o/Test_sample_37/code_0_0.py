# Initial matrix
matrix = [
    [0, 0, 0],
    [4, 32, 0],
    [0, 4, 0]
]

# Move sequence
sequence = 'uuu'

# Initial position of the number to move
row, col = 2, 1  # (3, 2) in 0-indexed

# Process each move in the sequence
for move in sequence:
    if move == 'u':
        if row > 0:
            # Check if the number can combine with the one above
            if matrix[row-1][col] == matrix[row][col]:
                matrix[row-1][col] *= 2
                matrix[row][col] = 0
            elif matrix[row-1][col] == 0:
                matrix[row-1][col] = matrix[row][col]
                matrix[row][col] = 0
            row -= 1

# Print the final matrix
print(matrix)