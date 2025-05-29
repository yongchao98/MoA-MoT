# Initial matrix
matrix = [
    [8, 0, 4],
    [32, 16, 16],
    [0, 0, 2]
]

# Move sequence
sequence = 'lll'

# Initial position of the number to move
row, col = 2, 2  # 0-indexed

# Process each move in the sequence
for move in sequence:
    if move == 'l' and col > 0:
        # Move left
        matrix[row][col - 1], matrix[row][col] = matrix[row][col], 0
        col -= 1

# Print the final matrix
print(matrix)